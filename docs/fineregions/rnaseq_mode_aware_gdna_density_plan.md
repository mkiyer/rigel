# RNA-Seq Mode-Aware gDNA Density Modeling Plan

Date: 2026-05-23
Status: planning document
Related docs: `fine_region_phase4_calibration.md`, `fractional_accumulator_phase3_plan.md`

## 1. Purpose

Phase 4 needs a gDNA estimator that solves the lab's most pressing use case
first: strand-specific hybrid capture RNA-seq. The broader architecture should
still support no-capture and unstranded libraries later, but the first useful
model should exploit the strongest orthogonal evidence available today:
strand-specific orientation.

The central separation remains:

1. The scanner and region accumulator collect mode-invariant evidence.
2. A library-mode policy decides which pieces of evidence are trustworthy.
3. A strand/depth model converts trustworthy evidence into expected gDNA mass,
   an upper gDNA bound, and a lower confident RNA bound.

For strand-specific hybrid capture, the first model can deconvolve observed
unspliced counts directly from strand balance. Capture-aware density and
exposure are still needed later, but they are not required to answer the first
region-level question: "given the observed sense/antisense split, how many RNA
fragments can we assert with high confidence?"

## 2. Library Modes And Priority

The two major axes are capture enrichment and strandedness:

| Mode | Difficulty | Main signal | Recommended status |
| --- | --- | --- | --- |
| Strand-specific, hybrid capture | Important and tractable | Capture-enriched regions plus RNA-vs-gDNA strand contrast | First implementation target |
| Strand-specific, no capture | Easiest | Same strand contrast plus cleaner background density | Reuse the strand model after capture target |
| Unstranded, no capture | Moderate | Intergenic and intron-only density with contained and boundary-side evidence | Later density-model target |
| Unstranded, hybrid capture | Hardest | Needs capture metadata or weak bounds; exon RNA and captured gDNA are not identifiable from counts alone | Conservative fallback only |

This table implies that Rigel should not have one monolithic "region density"
algorithm. It should have one evidence table and several estimator policies that
emit the same downstream products.

## 3. Core Architecture

### 3.1 Measurement layer: mode-invariant evidence

The C++ accumulator and Python payload views should not know whether a library
is capture or non-capture. They should continue to expose primitive
measurements:

- region geometry: `ref_id`, `start`, `end`, `signature`, adjacent signatures;
- unspliced contained counts by genomic strand;
- unspliced left-boundary and right-boundary counts by genomic strand;
- spliced contained and boundary counts by genomic strand for RNA diagnostics;
- contained effective length under the gDNA FL distribution;
- boundary-side effective exposure under the gDNA FL distribution;
- later: mappability-corrected and capture-corrected exposure fields.

Do not destructively combine contained and boundary evidence in this layer. A
combined count is a policy decision, not a storage decision.

An optional `EdgeView` can make boundary evidence easier to reason about, but
the logical output still needs side-specific mass because the left and right
sides of a junction can have different signatures and, under capture, different
targetability.

### 3.2 Library mode layer

Introduce an explicit library-mode object, even if the first implementation
only supports strand-specific hybrid capture:

```python
@dataclass(frozen=True, slots=True)
class LibraryMode:
    strandedness: Literal["unstranded", "strand_specific"]
    capture: Literal["none", "known", "inferred", "unknown"]
    capture_bed: Path | None = None
```

Strandedness is inferred from the existing global strand model trained on
spliced reads. Capture should eventually support both explicit user input and
conservative auto-detection, but the first capture-aware implementation should
prefer an explicit capture BED or panel definition.

### 3.3 Common downstream contract

Every mode should ultimately produce the same downstream table:

```python
@dataclass(frozen=True, slots=True)
class RegionGdnaEstimate:
    mean_count: np.ndarray       # expected gDNA fragments per region
    upper_count: np.ndarray      # high quantile, e.g. 95th percentile
    rna_lower_count: np.ndarray  # confident lower bound on RNA fragments
    exposure: np.ndarray         # mode-specific effective exposure, if known
    precision: np.ndarray        # prior strength / confidence proxy
    flags: np.ndarray            # fallback, capture, strand, low exposure, etc.
```

The per-locus EM prior aggregates `mean_count`, `upper_count`, and `precision`
over the locus footprint. It should not need mode-specific branching.

Important rule: a gDNA prior is uncertain evidence, not a hard physical count.
The precision field controls how strongly it enters EM.

## 4. Strand-Specific Hybrid Capture: First Model

### 4.1 Orientation convention

The existing `StrandModel` estimates:

```text
p = p_r1_sense = P(read 1 aligns in the transcript-sense direction)
```

The signed protocol value described in this plan is:

```text
SS = P(read 1 aligns in the transcript-antisense direction) = 1 - p
```

With this convention, `SS = 0` is R1-sense/R2-antisense, `SS = 0.5` is
unstranded, and `SS = 1` is R1-antisense/R2-sense. This signed value is not the
same as the current code's unsigned `StrandSummary.strand_specificity`, which is
`max(p, 1 - p)` and lives in `[0.5, 1.0]`.

The strand-deconvolution model needs both quantities:

```text
protocol_direction = sign(SS - 0.5)  # or equivalently sign(0.5 - p)
q = max(SS, 1 - SS)                  # RNA-major probability
```

Different stranded library preps rotate the informative strand in opposite
directions, so using only the unsigned specificity would be wrong.

For every non-ambiguous transcript-strand region (`TS_POS` or `TS_NEG`), first
compute transcript-relative counts with `sense_antisense_split(...)`. Then
rotate them into protocol-relative orientation:

```text
if p >= 0.5:
    RNA-major = transcript-sense
    RNA-minor = transcript-antisense
else:
    RNA-major = transcript-antisense
    RNA-minor = transcript-sense
```

After this rotation, RNA is expected on the RNA-major channel with probability
`q = max(SS, 1 - SS)`, regardless of whether the library is R1-sense or
R1-antisense. gDNA remains 50/50 in either orientation.

Regions with ambiguous transcript strand are excluded from this strand model.
Intergenic regions can still train the gDNA 50/50 strand variance using genomic
POS/NEG counts, but they cannot be used for transcript-relative RNA calls.

### 4.2 Generative model for one region

For a strand-informative region, define:

```text
k = observed RNA-major unspliced count
a = observed RNA-minor unspliced count
n = k + a
R = unknown RNA-derived unspliced count
D = unknown gDNA-derived unspliced count = n - R
q = max(SS, 1 - SS) = RNA-major probability for RNA fragments
```

The strand model is:

```text
DNA_major | D ~ BetaBinomial(D, mean=0.5, concentration=kappa_D)
RNA_major | R ~ Binomial(R, q)
k = DNA_major + RNA_major
```

The mean is:

```text
E[k | R] = 0.5 * (n - R) + q * R
         = 0.5 * n + (q - 0.5) * R
```

So the moment estimate is:

```text
R_hat = (k - 0.5 * n) / (q - 0.5)
```

clipped to `[0, n]`. For a nearly perfect stranded library (`q ~= 1`), this
reduces to:

```text
R_hat ~= k - a
D_hat ~= 2 * a
```

This is the intuitive result: the RNA-minor channel estimates half of the DNA,
and the RNA-major excess estimates RNA.

The variance is not modeled by strand fraction alone. Conditional on `R`:

```text
Var(DNA_major | D) = 0.25 * D * (D + kappa_D) / (1 + kappa_D)
Var(RNA_major | R) = R * q * (1 - q)
Var(k | R)         = Var(DNA_major | D) + Var(RNA_major | R)
```

This preserves the key dependence on total count. A 35/15 split with `n=50`
and an 800/200 split with `n=1000` do not carry the same certainty.

RNA can later be upgraded from binomial to beta-binomial if spliced-region
diagnostics show overdispersion. The first model should assume binomial RNA;
the protocol strand model is trained from many spliced reads and is usually
very precise.

### 4.3 Training the gDNA strand-balance overdispersion

gDNA has a fixed biological mean of 0.5. Only its overdispersion is learned.
For training region `i`, with counts `(x_i, n_i - x_i)` in any consistent
orientation, estimate `kappa_D` from the symmetric beta-binomial moment formula:

```text
S = sum_i (x_i - 0.5 * n_i)^2
B = 0.25 * sum_i n_i
A = 0.25 * sum_i n_i^2
kappa_D = (A - S) / (S - B)
```

with clipping to a finite range and a high-kappa fallback when `S <= B`.
This is the same statistical shape already used by the current strand-balance
diagnostic.

Training sources, in conservative order:

1. Intergenic regions: sparse in capture data, but high-purity gDNA. Use genomic
   POS/NEG counts because there is no transcript strand.
2. Intron-only regions with unique transcript strand: sparse but useful. Rotate
   into RNA-major/minor orientation when a transcript strand exists.
3. Intergenic-side and intronic-side boundary observations near exons: valuable
   in capture because they can be enriched while remaining gDNA-compatible.
4. Exons with no RNA-expression evidence: potentially the most useful captured
   training set, but must be selected conservatively.

Training regions must exclude ambiguous transcript-strand regions for any
transcript-relative operation. Ambiguous signatures need a density/imputation
path, not the strand deconvolution path.

Because capture changes how many fragments land in a region but not the 50/50
gDNA mean, captured DNA-rich regions are excellent for estimating strand
overdispersion. This is why the strand model can be implemented before the full
capture exposure model.

### 4.4 Selecting no-RNA exons for training

The statement "RNA-minor > RNA-major implies no RNA" is directionally right in
a highly stranded library, after protocol rotation. It should become a
statistical screen rather than a hard rule.

Candidate exon training regions must satisfy:

- unique transcript strand (`TS_POS` or `TS_NEG`), not ambiguous;
- low or absent spliced RNA evidence in the same region;
- no strong RNA-major excess in unspliced counts;
- enough total unspliced mass to be informative;
- optional: not overlapped by opposite-strand annotation or other ambiguity
  flags.

Initial practical rule:

1. Train `kappa_D` on the high-purity seeds above.
2. For each candidate exon, compute the posterior distribution over `R` using
   the mixture likelihood in Section 4.5.
3. Add the exon to the gDNA training set only if, for example,
   `P(R > max(1, 0.05 * n) | k, n) < 0.01`.
4. Refit `kappa_D` with the expanded training set.
5. Repeat once at most; avoid unstable self-training loops.

This includes the obvious `RNA-minor >= RNA-major` cases, but also permits
captured exons with a modest RNA-major excess when that excess is still
plausible under the learned gDNA beta-binomial.

Source-specific diagnostics should be kept: intergenic, intronic, boundary, and
no-RNA-exon training sets may have different apparent overdispersion. The
production model should prefer a conservative `kappa_D`, for example a pooled
estimate checked against the most overdispersed source stratum.

### 4.5 Posterior deconvolution and confident RNA lower bound

For one region, compute:

```text
P(R = r | k, n) proportional to P(k | R = r, n, q, kappa_D) * P(R = r)
```

where:

```text
P(k | R = r, n, q, kappa_D)
  = convolution(
        BetaBinomial(D = n - r, mean=0.5, kappa_D),
        Binomial(R = r, q)
    )
```

The first prior should be weak and conservative, such as uniform over `r` from
0 to `n`. Later, spliced counts can inform the prior that RNA expression is more
likely, but the lower-bound call should not depend on an aggressive expression
prior.

Outputs for confidence level `alpha`:

```text
RNA_lower_alpha = posterior_quantile(1 - alpha) of R
gDNA_upper_alpha = n - RNA_lower_alpha
RNA_mean = posterior_mean(R)
gDNA_mean = n - RNA_mean
```

The framing becomes:

```text
minimum RNA molecules asserted with alpha confidence = RNA_lower_alpha
maximum gDNA molecules still plausible at alpha confidence = gDNA_upper_alpha
```

This is better than regressing `log(antisense)` on `log(sense)` because it uses
the known component means, the measured protocol strand specificity, the
region's total count, and the learned DNA overdispersion directly.

Implementation note: region counts are fractional because the accumulator
partitions fragments across contained and boundary compartments. The initial
implementation can evaluate the posterior on an integer grid using rounded
effective counts for moderate/high counts and use a moment-matched normal
likelihood for low or highly fractional counts. Low-count regions should carry
`FLAG_LOW_EXPOSURE` and naturally yield weak lower RNA bounds.

### 4.6 Conservative tail approximation

The posterior grid is the principled path, but a fast conservative approximation
is useful for vectorized screening.

For each candidate RNA count `r`, compute the mean and variance in Section 4.2
and evaluate a normal or saddlepoint tail for the observed `k`. The smallest
`r` whose tail is compatible gives a lower confidence bound on RNA. This is a
screening approximation; final per-region outputs should use the posterior
where feasible.

In the near-perfect strand case, a helpful intuition remains:

```text
gDNA is given credit for a high-tail RNA-major imbalance;
only RNA-major counts above that DNA upper tail are called confident RNA.
```

### 4.7 Region outputs

The first strand-specific hybrid capture estimator should emit per-region
arrays:

- `n_unspliced_total`
- `rna_mean_count`
- `rna_lower_count_95` and optionally 99/99.9
- `gdna_mean_count`
- `gdna_upper_count_95` and optionally 99/99.9
- `strand_deconv_precision`
- `source_flags`

For EM prior assembly, the direct deconvolved gDNA count is already useful.
Capture-aware density is still needed for extrapolating to regions with little
or no observed mass, but it does not have to block the first strand model.

### 4.8 Limits

This model does not solve:

- transcript-strand-ambiguous regions;
- regions where both plus- and minus-strand annotations overlap;
- unstranded hybrid capture;
- expression generated by unannotated opposite-strand transcripts;
- density extrapolation into completely unobserved regions.

Those cases need density-based estimation, capture-aware exposure, neighbor
imputation, or conservative fallback flags.

## 5. Boundary Evidence In The Strand Model

For strand-specific hybrid capture, boundary evidence should be source-aware.

Useful sources:

- intronic or intergenic side of exon-intron/exon-intergenic boundaries:
  gDNA-compatible and often capture-enriched;
- exon side of a boundary: usable only after the no-RNA or posterior
  deconvolution screen;
- contained exon mass: usable for deconvolution, not automatically for training.

The model can run on any source that provides `(RNA-major, RNA-minor, total)`
counts and a non-ambiguous transcript strand. Boundary-side storage should not
be collapsed before the evidence policy decides which side is being modeled.

## 6. Capture Exposure And Density

The strand deconvolution model estimates observed gDNA count. A capture-aware
density model is still needed for:

- predicting gDNA in low-count or unobserved regions;
- constructing regional exposure weights;
- comparing across targeted and off-target regions;
- building priors when local observed evidence is weak.

With a capture BED or BED12 panel, build a capture targetability field and use
it in the exposure denominator:

```text
mu_R = rho_background * E_background_R + rho_capture * E_capture_R
```

or equivalently:

```text
mu_R = rho_g * (E_background_R + eta_capture * E_capture_R)
```

where:

- `E_background_R` is FL-corrected opportunity outside target influence;
- `E_capture_R` is FL-corrected opportunity inside target influence;
- `eta_capture` is the enrichment ratio, learned or supplied;
- both contained and boundary-side exposures need capture-aware versions.

The capture exposure model should be trained from deconvolved gDNA counts, not
raw unspliced exon counts.

## 7. No-Capture Estimator

For no-capture libraries, the earlier density plan remains useful but is no
longer the first implementation target.

For each gDNA-compatible region side:

```text
N = unspliced_contained + compatible_unspliced_boundary_left + compatible_unspliced_boundary_right
E = contained_effective_length + compatible_boundary_left_exposure + compatible_boundary_right_exposure
```

Then fit a negative-binomial density model with a quantile interface:

```text
N_R ~ NB(mu_R, kappa_R)
log(mu_R) = beta_0 + beta_1 * log(E_R)
log(kappa_R) = gamma_0 + gamma_1 * log(E_R)
```

The current pooled density estimator is the special case `beta_1 = 1` and
constant `kappa`. This model remains a good later target for no-capture and for
ambiguous regions that cannot use strand deconvolution.

## 8. Implementation Roadmap

### M1. Protocol-rotated strand observation builder

Build a vectorized helper that returns RNA-major/RNA-minor counts for each
non-ambiguous region and each source compartment:

- contained unspliced;
- left-boundary unspliced;
- right-boundary unspliced;
- contained and boundary spliced diagnostics.

The helper must use the signed `p_r1_sense`, not only `SS`, so R1-sense and
R1-antisense preps are handled correctly.

### M2. gDNA strand-balance training set

Implement high-purity seed extraction and the symmetric beta-binomial `kappa_D`
fit. Keep source-specific diagnostics and exclude ambiguous transcript-strand
regions.

### M3. No-RNA exon selection

Use the seed `kappa_D` and spliced evidence to select captured exons that are
consistent with no RNA expression. Refit `kappa_D` once with these additional
training regions.

### M4. Posterior strand deconvolution

Implement per-region posterior deconvolution for `R`, returning RNA lower bounds
and gDNA upper bounds at 95, 99, and optionally 99.9 percent confidence.

### M5. RegionGdnaEstimate integration

Store strand-deconvolved gDNA mean/upper counts and precision in the common
`RegionGdnaEstimate` table. Use these counts directly for the first EM prior
restoration path.

### M6. Capture exposure model

Add optional capture BED/BED12 ingestion and compute capture-corrected exposure.
Train density/exposure models from deconvolved gDNA counts, not raw counts.

### M7. Density models for fallback cases

Implement no-capture NB density modeling and ambiguous-region fallback. Use
these paths for unstranded data and transcript-strand-ambiguous regions.

## 9. Validation Plan

Each milestone needs simulations that isolate the intended signal:

- R1-sense and R1-antisense strand-specific libraries: verify protocol rotation
  produces identical deconvolution after rotation;
- strand-specific hybrid capture, gDNA-only targeted exons: recover 50/50 mean
  and estimate realistic `kappa_D`;
- captured expressed exons: RNA lower bound should be positive only when the
  RNA-major excess exceeds the learned DNA upper tail;
- low-count exons: posterior intervals should widen and lower RNA bounds should
  often be zero;
- exons with spliced evidence: should not enter no-RNA training, but should be
  deconvolved for output;
- ambiguous signatures: should be excluded from strand deconvolution and routed
  to density/imputation fallback;
- capture BED simulations: deconvolved gDNA counts should train capture exposure
  without raw RNA contamination.

Acceptance metrics should include calibration of RNA lower bounds, coverage of
true gDNA counts by `gdna_upper_count_alpha`, source-specific `kappa_D`
stability, and downstream mRNA/nRNA/gDNA assignment error in EM.

## 10. Near-Term Recommendation

Start with the strand-specific hybrid capture model, not the no-capture density
model. In practice:

1. Build protocol-rotated RNA-major/RNA-minor observations.
2. Train a conservative symmetric beta-binomial gDNA strand-balance model.
3. Select no-RNA captured exons only through a posterior screen.
4. Deconvolve each non-ambiguous region into RNA lower bound and gDNA upper
   bound.
5. Store direct deconvolved gDNA counts for EM prior restoration.
6. Add capture-aware exposure afterward for density extrapolation and low-count
   fallback.

This path attacks the lab's pressing capture use case first while preserving a
general architecture for no-capture, unstranded, and ambiguous-region models.