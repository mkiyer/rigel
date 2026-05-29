# Hybrid-Capture Failure Research Memo

Date: 2026-05-25
Status: research / design memo
Related benchmark report: `docs/benchmarks/hyb_capture_500kb_calibration_eval_2026-05-25.md`
Related design proposal: `docs/calibration/hybrid_capture_exposure_model_2026-05-18.md`
Diagnostic script: `scripts/debug/hybrid_capture_failure_audit.py`

## 1. Short Answer

The hybrid-capture failure is not primarily a scanner, accumulator, or EM plumbing
failure. The foundation is mostly the one we wanted: fine regions, fractional
contained/boundary counts, FL-aware effective opportunities, a Gamma/negative-binomial
density model, and strand-specific regional deconvolution all exist.

The problem is that the current density model is still effectively an off-target
background model. It treats the low intron/intergenic density seen under capture as a
strong prior, then asks local boundary flux to update that prior. In capture libraries
that prior becomes both low and extremely inert, so boundary evidence cannot move the
regional exposure enough. The resulting `RegionExposure.A_r` only reaches about 1.9 in
the synthetic capture-high run, while the simulator creates a several-hundred-fold
captured/uncaptured transcript imbalance.

Separately, transcript quantification under capture is missing an RNA exposure model.
Even with no gDNA, uncaptured transcripts lose essentially all observed fragments and
Rigel converts that missing exposure into near-zero TPM. That all-transcript TPM is not
identifiable from the BAM alone unless Rigel has, or learns with strong assumptions, a
capture target/exposure model.

## 2. Evidence From The Current Suite

`hyb_capture_500kb` has 8 conditions: capture on/off, ss=0.99/0.50, gDNA none/high.
The compact audit script reproduces the key deltas from existing outputs.

| Metric | gDNA high ss=0.99 capture off | gDNA high ss=0.99 capture on | Interpretation |
|---|---:|---:|---|
| `rho_ref` | 0.2002 | 0.0029 | off-target anchor collapses 69x |
| intron contained FL-pool mass | 27,588 | 394 | capture removes intronic contained gDNA |
| intron boundary FL-pool mass | 1,422 | 4,126 | boundary evidence is enriched and should be useful |
| exon contained FL-pool mass | 84,646 | 158,659 | exon/capture region receives the mass |
| exon boundary FL-pool mass | 1,719 | 12,422 | captured exon boundaries light up |
| gDNA FL mean | 152.1 | 168.3 | observed gDNA is capture-conditioned / length-biased |
| `RegionExposure.A_max` | 1.07 | 1.91 | exposure learner sees only a 2x effect, not capture-scale enrichment |
| sum gDNA prior count | 35,820 | 17,121 | prior mass goes down while observed gDNA pileup goes up |

At the locus level, capture-on gDNA counts increase, but priors shrink:

| Example locus | `gdna_off` | `gdna_on` | `gdna_prior_off` | `gdna_prior_on` | prior fold |
|---:|---:|---:|---:|---:|---:|
| 1 | 6,955 | 16,243 | 0.092 | 0.034 | 0.37 |
| 8 | 6,669 | 15,595 | 0.953 | 0.134 | 0.14 |
| 4 | 3,842 | 10,667 | 0.780 | 0.112 | 0.14 |
| 9 | 594 | 4,634 | 0.137 | 0.031 | 0.23 |

The transcript exposure failure is independent of gDNA. Expressed transcripts with any
probe coverage carry 646k truth TPM; uncaptured expressed transcripts carry 353k truth
TPM. Under capture-on/no-gDNA, Rigel reports approximately 999k TPM for captured
transcripts and only 670-860 TPM for uncaptured transcripts. This proves the abundance
collapse is a capture-denominator problem, not a contamination problem.

## 3. How The Code Matches The Original Concept

The code matches several major parts of the roadmap.

1. The C++ accumulator preserves the measurement layer we wanted. It records 12-channel
   fractional counts per fine region: contained/boundary-left/boundary-right x
   spliced/unspliced x pos/neg. It enforces one fragment equals one unit of fractional
   mass.
2. `DensityObservation` keeps contained and boundary evidence separate and computes
   separate FL-aware contained and boundary opportunities.
3. `density_model.py` implements the intended conjugate skeleton:
   `alpha_post = alpha0 + B`, `beta_post = beta0 + L_b`, then a negative-binomial
   predictive distribution over contained gDNA.
4. `strand_deconv.py` implements the strand-specific question from the TODO: given
   sense/antisense counts, infer a posterior over RNA and gDNA using a binomial RNA
   component and beta-binomial DNA strand-balance component.
5. `integration.py` fuses density and strand evidence into bounded regional gDNA
   posteriors.
6. `prior.py` projects regional fused gDNA evidence to MultiLocus prior counts and
   computes a gDNA effective-length denominator.

So the project did not go astray at the infrastructure level. The measurement and
posterior framework are pointed in the right direction.

## 4. Where We Went Astray

### 4.1 Boundary evidence is present but mathematically neutered under capture

The intended idea was: off-target regions shrink to the global background, while
on-target regions with large boundary flux overwhelm the prior and snap to local
exposure.

In the current implementation, the prior precision cap produces `alpha ~= 400` for
successful anchor fits. In capture-on/high-gDNA, the intron mean density is tiny:

```text
rho ~= 0.0029
alpha0 ~= 400
beta0 = alpha0 / rho ~= 136,000 bp
```

The local update is then:

```text
rho_post = (alpha0 + B) / (beta0 + L_b)
```

A typical boundary opportunity `L_b` is tiny relative to `beta0`, so the posterior is
prior-dominated even when boundary counts are biologically meaningful. The summary says
`n_prior_dominated = 125 / 125` in the failing capture-high stranded condition.

This is the main mathematical mismatch. The code implements the Gamma update, but the
prior concentration is not appropriate for learning capture exposure from sparse local
boundary evidence.

### 4.2 The model has one background density, not a capture exposure model

Hybrid capture is bimodal or multi-modal: off-target introns/intergenic regions are
depleted, while targeted exon-adjacent regions are enriched. The current model has
family priors (`INTERGENIC`, `INTRON`, `ALL`) plus per-region updates, but it does not
fit an explicit on-target/off-target exposure mixture. It therefore interprets the
low intron density as low library-wide gDNA availability rather than low off-target
exposure.

The synthetic data is screaming this distinction:

```text
intron contained: 27,588 -> 394
intron boundary:   1,422 -> 4,126
exon boundary:     1,719 -> 12,422
```

Contained introns and exon boundaries move in opposite directions. The current density
model compresses that into a single low reference density and a weak local A-ratio.

### 4.3 Exposure is only a gDNA denominator, and only at locus resolution

The current live quant path uses `RegionExposure` to compute a bp-weighted mean
exposure over each MultiLocus and multiplies the gDNA effective length by that mean.
That is incomplete for capture.

The generative model in the earlier design requires both sides:

```text
gDNA score(f) += log W(f)
gDNA denominator = integral W(s, ell) over valid starts
```

A locus-wide mean is especially weak when a large locus has small captured islands.
The high-exposure fragments need local `log W(f)` support, not only a diluted locus
average. The current data show locus `gdna_eff_len_weight_ratio` only around 1.2-1.4
in capture-on, far below the target/off-target contrast.

### 4.4 RNA exposure is absent

This is the catastrophic transcript-abundance failure. The ordinary RNA denominator is
still transcript effective length, not capture-accessible transcript effective length.
When capture removes uncaptured transcripts from the observed library, Rigel interprets
that as low expression.

Without a probe BED or a learned transcript exposure model, all-transcript TPM for
untargeted or weakly targeted transcripts is not identifiable from the BAM alone. Rigel
can still produce useful target-limited quantification and clear low-exposure flags, but
it should not silently report uncaptured transcripts as true near-zero expression.

### 4.5 Strand evidence is useful but cannot rescue the whole capture problem

Strand-specific deconvolution is the best near-term lever for ss=0.99 capture data. It
can estimate how much unspliced mass is DNA-like versus RNA-like in strand-informative
regions. But today it is fused with a sharp, wrong density prior, and it does not solve:

- unstranded capture (`ss=0.50`);
- ambiguous transcript-strand regions;
- transcript-level RNA exposure;
- the gDNA local scoring denominator.

It should become the primary gDNA numerator signal for strand-specific capture, while
boundary exposure supplies the denominator and handles unstranded/ambiguous cases.

### 4.6 Fragment-length behavior is being mixed conceptually

The observed capture-on gDNA FL mean is 168 bp vs simulated pre-capture mean 150 bp.
That is not simply an estimation bug: capture selection makes longer gDNA fragments more
likely to overlap a probe. For scoring observed captured fragments, the conditional
observed gDNA FL may actually be the right distribution. For estimating latent library
contamination or effective opportunity, we need a selection-aware denominator.

So the fix is not merely "force gDNA FL back to 150." The cleaner approach is to keep
separate concepts:

```text
latent molecule FL distribution
observed capture-conditioned FL distribution
capture-weighted effective opportunity
```

## 5. Solution Set To Consider

### Solution 1: Add a capture-mode detector immediately

This is low risk and high diagnostic value. Capture is obvious in the current outputs:

- intron/intergenic contained density collapses;
- exon boundary and exon contained mass increase;
- `tail_probability` and `expected_tail_count` rise;
- `A_max` remains low despite boundary enrichment;
- gDNA FL shifts upward.

Emit a `capture_mode: uniform | suspected_capture | capture` diagnostic in
`summary.json`. Use it to choose safer calibration paths and separate acceptance
thresholds.

A simple first detector can be rule-based:

```text
capture_suspected if
  exon_boundary_mass / max(intron_boundary_mass, 1) is high
  and intron_contained_density << exon_boundary_density proxy
  and density_evidence.n_prior_dominated is large
```

Later this can become a small likelihood-ratio classifier.

### Solution 2: Replace the current local exposure update with a boundary exposure model

Keep the conjugate Gamma/NB core, but change what it represents.

Use contained intergenic/intronic regions to estimate `rho_off`, the off-target gDNA
background density. Then fit boundary rates separately:

```text
B_e ~ NegBin / Gamma-Poisson(lambda_e * L_b)
lambda_e = rho_off * A_e
A_e drawn from a small exposure mixture: off-target, mid, on-target
```

The output is a posterior over local exposure:

```text
E[A_e | B_e]
P(on-target | B_e)
```

Then predict contained gDNA using the mixture posterior:

```text
D_contained ~ mixture_z NegativeBinomial(alpha_z + B, beta_z + L_b, L_c)
```

Key difference from the current code: the local exposure prior should have modest
concentration. It should not use `alpha ~= 400` from the off-target contained-density
fit. The background density can be well-estimated globally while local exposure remains
uncertain and movable.

For non-capture libraries, the fitted mixture should collapse to one class and return
`A_r ~= 1`.

### Solution 3: For strand-specific capture, let strand deconvolution own the numerator

In `ss ~= 1` capture data, strand is orthogonal evidence. Use it more directly:

1. Estimate `kappa_D` from pure gDNA seeds plus self-trained no-RNA exons.
2. For strand-informative regions, compute the posterior over `D_r` from strand counts.
3. Use density/boundary exposure mainly as a denominator and uncertainty prior, not as a
   sharp low-count veto.
4. If density and strand disagree strongly in capture mode, emit a tension flag and
   prefer the strand-derived gDNA upper/mean for the EM prior numerator.

This is the fastest path to improving the most important mode: strand-specific hybrid
capture.

### Solution 4: Restore local exposure in gDNA scoring, not only locus effective length

The gDNA model should use the same exposure field in the numerator likelihood and the
denominator:

```text
log_lik_gdna(f) += log W(f)
L_gdna(locus) = sum_ell h_G(ell) * integral valid_start W(s, ell)
```

The current locus bp-weighted mean exposure is too diluted. A captured exon inside a
large locus needs the fragment-level `log W(f)` boost at its genomic midpoint or footprint.

Implementation shape:

- introduce a `CaptureExposureModel` / `RegionExposure` with piecewise-constant intervals;
- add `log_weight_for_fragment_footprint` or `log_weights_for_midpoints`;
- apply this before partitioning EM data;
- replace the current mean-weighted denominator with an FL-convolved weighted integral.

This aligns with the existing design proposal and the earlier regional-exposure memory,
but the live path appears to have lost the per-unit scoring application.

### Solution 5: Add RNA exposure semantics before trusting TPM under capture

There are two honest modes.

Known-target mode:

- accept a probe BED/transcript-probe file;
- build transcript capture effective lengths;
- use capture-weighted RNA denominators for mRNA;
- use capture-weighted genomic denominators for gDNA/nRNA.

No-target mode:

- infer targeted regions from boundary/exposure evidence;
- report target-limited or exposure-conditioned quantification;
- flag transcripts with insufficient exposure as `low_exposure` / `not_quantifiable`;
- avoid interpreting no observed fragments as true zero expression.

For fully untargeted transcripts, all-transcript TPM is not identifiable from a single
capture BAM. The model can learn that they are depleted; it cannot recover their true
abundance without external target geometry, a strong abundance prior, or cross-sample
information.

### Solution 6: Treat gDNA priors as uncertain evidence

Current `gdna_prior_count_em` is a hard physical pseudocount. The TODO already flags
this. We should carry precision/variance into EM or at least into prior scaling.

Practical first options:

- use posterior mean plus a concentration scalar rather than raw expected count;
- cap prior strength when density and strand disagree;
- in capture mode, use `upper_count` as a conservative gDNA-protective prior only when
  no-gDNA controls do not inflate false gDNA;
- expose `prior_precision` in `loci.feather` for debugging.

Longer-term, integrate the regional gDNA posterior as a hierarchical prior instead of a
single count.

### Solution 7: Fix evaluation gaps

The current synthetic suite is still useful, but it lacks robust intergenic anchor
coverage. Add:

- realistic intergenic padding and mappability masks;
- nRNA + capture combinations;
- known-probe oracle runs;
- strand-only, density-only, and exposure-only ablations;
- gDNA leak stratified by captured exon, uncaptured exon, intron, intergenic, and
  transcript-strand ambiguity;
- persisted per-region density/strand/exposure tables for debugging.

## 6. Recommended Implementation Path

### Phase A: Stop the bleeding in strand-specific capture

1. Add capture detector diagnostics.
2. In capture-detected, highly stranded libraries, reduce density-prior dominance for
   exonic/non-anchor regions and use strand-derived gDNA evidence more directly.
3. Add summary fields that show density-vs-strand tension and prior-dominated regions.
4. Validate on the 8-condition suite. Target: gDNA->RNA leak in ss=0.99 capture-high
   drops substantially without creating gDNA in capture-on/no-gDNA.

This does not solve transcript TPM yet, but it should improve classification.

### Phase B: Boundary-driven exposure model

1. Build an edge/boundary table from existing left/right region counts.
2. Fit an off-target background from contained anchors.
3. Fit a two- or three-component boundary exposure mixture.
4. Produce `A_r`, `P_on`, and uncertainty per region.
5. Fuse density/strand using this exposure posterior.
6. Restore local `log W(f)` scoring and weighted gDNA effective length.

Target: `A_max` and on/off contrast become capture-scale, per-locus priors at captured
loci rise instead of collapsing, and gDNA->RNA leak approaches capture-off levels.

### Phase C: RNA/transcript exposure

1. Add optional known-probe path first because it is identifiable and testable.
2. Compute capture-weighted mRNA effective lengths.
3. Add no-target inferred exposure flags and target-limited quant semantics.
4. Only then evaluate all-transcript TPM under capture.

Target: capture-on/no-gDNA no longer reports uncaptured transcripts as confidently zero;
known-probe runs recover captured/uncaptured aggregate TPM much more honestly.

## 7. Design Decisions I Recommend

1. Keep the current fine-region/fractional accumulator. It is the right substrate.
2. Split `DensityEvidence` from `ExposureModel`. Density is contamination rate;
   exposure is assay opportunity. They are currently entangled.
3. Do not use the off-target Gamma precision (`alpha ~= 400`) as the local exposure
   prior precision. This is the immediate reason boundary evidence cannot move.
4. Prioritize strand-specific hybrid capture first. It is both the important use case
   and the one with the strongest orthogonal evidence.
5. Be explicit about identifiability: no-BED capture can learn target/off-target
   structure well enough for gDNA calibration, but it cannot fully recover expression
   for transcripts with no exposure.
6. Treat capture-conditioned gDNA FL as an observed-library model, not simply as a bad
   estimate. Correct denominators and selection effects rather than forcing the FL mean
   to the pre-capture simulator parameter.

## 8. Bottom Line

We did not take a wrong architectural road. We built most of the hard measurement
infrastructure. The current failure comes from using that infrastructure with a density
model that still assumes one off-target background can be locally updated into capture
exposure. Under hybrid capture, the off-target prior becomes low, sharp, and wrong for
exonic probe regions.

The fix is to promote exposure to a first-class latent field: learned from boundary and
strand-compatible gDNA evidence, used consistently in gDNA scoring and denominators, and
then extended to RNA effective lengths or explicit low-exposure transcript semantics.
