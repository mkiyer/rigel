# Hybrid-Capture Exposure Model

Status: design proposal, not implemented

Date: 2026-05-18

## Problem

Rigel's current gDNA component assumes that genomic opportunity is close to
uniform after mappability and fragment-length geometry. That is reasonable for
polyA/ribo-minus style libraries where contaminating gDNA fragments sample the
mappable genome broadly. It is wrong for hybrid-capture libraries.

In capture data, the sampling process is inhomogeneous: on-target bases have
high fragment density and off-target bases can be depleted by orders of
magnitude. A `MultiLocus` can span 100 kb while only 10 kb are actually
captured. If the gDNA EM component is normalized by the whole 100 kb footprint,
Rigel divides the gDNA mass by a denominator dominated by off-target introns
that do not generate fragments. That artificially lowers the gDNA posterior for
on-target unspliced fragments and lets mRNA/nRNA components win.

The fix should not be a hard-coded capture-only mode. Rigel needs a single
exposure framework that behaves like the current uniform model for polyA and
ribo-minus data, but learns or accepts a nonuniform target/off-target profile
for hybrid capture.

Terminology note: nRNA is unspliced by definition. Spliced alignments are mature
RNA evidence or alignment artifacts, not nRNA evidence. Any capture-aware model
should preserve that invariant.

## Generative Model

Introduce a genomic exposure field `A(x)` describing relative sampling
opportunity at genomic position `x`. For a fragment of length `ell` starting at
`s`, define a fragment exposure weight:

```text
W(s, ell) = mean_or_integral_{x in [s, s + ell)} A(x)
```

For gDNA in a `MultiLocus` `M`, the likelihood becomes an inhomogeneous genomic
fragment process:

```text
P(f | gDNA_M) ∝ h_G(ell_f) q_G(f) W(f) / Ltilde_gDNA(M)

Ltilde_gDNA(M) = sum_ell h_G(ell) sum_s W(s, ell) I([s, s+ell) overlaps M)
```

`h_G` is the gDNA fragment-length PMF and `q_G(f)` contains non-exposure terms
such as mismatch, strand symmetry, and splice/gDNA eligibility. This is the
weighted version of the current `_exposure.gdna_eff_len_for_loci()` contract.

The exposure scale is not identifiable separately from the global gDNA density:
`A(x) * c` and `rho_gDNA / c` produce the same expected counts. Therefore only
the shape matters. Normalize `A` by a robust convention, for example median
off-target exposure = 1 or 95th-percentile exposure = 1. The same normalized
field must be used in both terms:

1. add `log W(f)` to the gDNA fragment score;
2. replace unweighted `gdna_eff_len` with weighted `Ltilde_gDNA(M)`.

If this is done consistently, uniform libraries naturally reduce to today's
model with `A(x) ≈ 1` everywhere.

## Why This Fixes The Current Failure

Suppose a capture locus spans 100 kb, with 10 kb on-target and 90 kb off-target.
If off-target exposure is 1/100 of on-target exposure, the correct gDNA
opportunity is approximately:

```text
10 kb * 1.0 + 90 kb * 0.01 = 10.9 kb
```

The current unweighted denominator is approximately 100 kb. That penalizes gDNA
by:

```text
log(100 / 10.9) ≈ 2.2 nats
```

before considering the fragment's own local target weight. This is large enough
to flip ambiguous unspliced fragments to RNA when the mRNA/nRNA component has a
large abundance/effective-length advantage. In larger multi-loci the penalty can
be substantially worse.

## Two Inputs, One Model

Rigel should support both supervised and unsupervised exposure construction.

### Optional BED Path

If a target BED is available, build `A(x)` from it:

```text
A(x) = a_on     inside target intervals and probe shoulders
A(x) = a_off    outside target intervals
```

The on/off ratio should still be estimated from the BAM during calibration,
not hard-coded. BED supplies the geometry; calibration estimates the fold
change and smoothing width. This avoids the common failure where a provided BED
is close but not exactly the assay-specific capture footprint.

### No-BED Path

When no BED is available, infer `A(x)` from the alignment itself using a robust
empirical-Bayes exposure model. This does not need to recover the exact bait
design. It only needs to learn enough of the target/off-target contrast to avoid
normalizing on-target gDNA by off-target intronic span.

## Unsupervised Exposure Learning

### 1. Build Exposure Tiles

Create a genome partition for exposure learning. This can reuse `index.region_df`
and add a coarser fixed-window layer.

Recommended first representation:

```text
tile = intersection(region_df segment, fixed windows of 500 bp to 2 kb)
```

Record for each tile:

- genomic interval and region type: exon, intron, intergenic;
- mappability or valid-start exposure if available;
- mature-RNA evidence: annotated spliced counts overlapping the tile;
- gDNA-proxy evidence: intergenic-contained, intron-contained, exon-intron
  boundary, opposite-strand unspliced in stranded libraries, and high gDNA
  posterior fragments from a preliminary pass;
- ambiguous unspliced counts, retained for later but not trusted initially.

The initial tile model should use conservative gDNA-proxy evidence. Ambiguous
exon-contained unspliced fragments are exactly the failure mode and should not
be allowed to define the target map on the first pass.

### 2. Fit A Small Mixture Of Exposure States

Model gDNA-proxy counts as overdispersed counts with tile-specific exposure
state:

```text
y_i ~ NegBin(mu_i, kappa)
mu_i = rho_{type(i)} * A_i * E_i
log A_i in {a_off, a_mid, a_on}       # 2-3 states are enough initially
```

`E_i` is the FL-aware valid-start exposure for the tile. `rho_type` preserves
the existing intergenic/intron/boundary density channels. The latent state
captures assay enrichment independent of region type.

Fit this with empirical Bayes:

1. initialize `A_i = 1`;
2. estimate `rho_type` from robust central quantiles of gDNA-proxy density;
3. fit a mixture to residual log density `log((y_i + c) / (rho_type E_i + c))`;
4. update posterior `E[A_i | y_i, features]`;
5. smooth/posterior-regularize neighboring tiles;
6. repeat once after a preliminary EM pass supplies better gDNA posterior
   counts.

For non-capture libraries, the fitted mixture should collapse: the estimated
on/off log contrast will be near zero or the posterior mass will concentrate in
one state. In that case Rigel should explicitly mark the exposure model as
`uniform` and use the existing unweighted path.

### Identifiability Boundary

Keep these quantities conceptually separate:

```text
A(x)       assay exposure / capture accessibility
rho_type   global gDNA contamination density after exposure correction
theta      locus EM abundance/mass parameter
```

The no-BED learner must not become an arbitrary local density smoother. If every
locus can get its own free `A`, then `A`, `rho`, and `theta` are not separately
identified and the model can explain away errors instead of fixing them. The
exposure field should therefore be low-complexity and assay-like: a small number
of states, shared genome-wide fold changes, weak spatial regularization, and
strong shrinkage to uniform when evidence is limited.

### 3. Use RNA Evidence As A Confounder, Not As Target Evidence

High mature-RNA expression can mimic capture enrichment if total exon coverage
is used naively. The no-BED learner should therefore treat spliced mature-RNA
coverage as a confounder:

```text
ambiguous_exon_count_i ≈ gDNA_i + mature_RNA_i + nRNA_i
```

For exposure learning:

- do not use annotated spliced mature-RNA counts as positive evidence that a
  region is targeted;
- use mature-RNA counts to downweight ambiguous exon-only counts during the
  first exposure fit;
- allow boundary-crossing and intron/intergenic gDNA-proxy counts to anchor the
  exposure state near expressed exons;
- after one EM pass, use posterior expected gDNA counts rather than total counts.

This is conservative. It may miss some targets, but it should avoid the more
dangerous error of interpreting expression as gDNA capture exposure.

### Avoid Circular Self-Rescue

The same ambiguous unspliced fragments should not be allowed to both define a
tile as on-target and then benefit from the reduced denominator caused by that
definition. Use at least one of these safeguards:

- learn `A(x)` from conservative gDNA-proxy categories before using ambiguous
  exon-contained fragments;
- cross-fit tiles, for example learn exposure on odd fragments and score even
  fragments, then swap;
- leave one locus or one genomic block out when estimating the exposure state
  used for that locus's gDNA denominator;
- cap the contribution of posterior-gDNA ambiguous fragments during refinement.

This is especially important in long nRNA loci, where the current failure mode
already has a self-seeding character.

### 4. Spatial Regularization

Capture exposure is locally structured but not necessarily smooth at exon
boundaries. Use weak regularization, not broad smoothing:

- merge adjacent tiles with similar posterior exposure;
- allow sharp exon/intron transitions;
- optionally add shoulders around high-exposure exons, learned from boundary
  events and fragment length;
- cap extreme on/off ratios unless many independent tiles support them.

The output should be a compact piecewise-constant `CaptureExposureModel`, not a
per-base array.

## Projecting The Exposure Field

### gDNA Effective Length

Replace `_exposure.gdna_eff_len_for_loci()` with a weighted sibling:

```text
weighted_gdna_eff_len_for_loci(loci, ref_lengths, fl, exposure_model)
```

For each fragment length `ell`, expand every locus interval into the valid-start
window `[start - ell + 1, end)`, clip to contig-valid starts, intersect with the
piecewise exposure intervals, merge overlapping windows, and sum weighted
lengths. The unweighted function is the special case `A(x)=1`.

For performance, precompute interval-prefix sums of `A(x)` per reference so a
weighted window integral is `O(log n)` in the number of exposure segments. The
current loop over positive fragment lengths can stay; only the `n_valid` count
changes from raw length to weighted length.

### gDNA Fragment Score

The scorer needs the same field:

```text
score_gDNA(f) += log W(f)
```

For a paired fragment, `W(f)` should be the average exposure over the inferred
genomic footprint, or the exposure integrated over its aligned blocks if the
insert is uncertain. A floor such as `W_min = 1e-3` prevents impossible calls in
off-target regions and keeps log likelihoods finite.

### gDNA Priors

`expected_gdna_count_global()` currently projects global densities onto
unweighted exon/intron/intergenic exposure. It should project onto weighted
exposure instead:

```text
eta_g(locus) = rho_ig * L_ig_weighted
             + rho_in * L_in_weighted
             + rho_boundary * (B_cross_weighted + L_exon_weighted)
```

The key invariant is unchanged: prior mass and EM effective length must describe
the same physical opportunity space.

### Mature RNA And nRNA

Hybrid capture affects RNA molecules too. The minimal fix can start with gDNA
because the observed failure is gDNA under-normalization, but the principled
model should eventually project the same exposure field onto RNA components:

- mature RNA: exposure-weighted effective length over spliced exonic
  coordinates;
- nRNA: exposure-weighted effective length over unspliced genomic transcript
  spans, with only unspliced alignments eligible;
- spliced alignments should never be treated as nRNA evidence.

This keeps transcript quantification sane for capture panels and avoids moving
the bias from gDNA into RNA effective lengths.

## Calibration Flow

Recommended first implementation sequence:

1. Run the current scan and SRD calibration to obtain fragment-length models,
   strand model, region counts, and preliminary assignments.
2. Build exposure tiles and conservative gDNA-proxy counts.
3. Fit `CaptureExposureModel` and classify the library as `uniform` or
   `inhomogeneous` based on fitted on/off contrast and held-out likelihood gain.
4. Recompute global gDNA densities using weighted denominators.
5. Recompute per-locus `gdna_prior_count` and weighted `gdna_eff_len`.
6. Re-score gDNA candidates with `log W(f)` and rerun locus EM.
7. Optionally run one refinement iteration using posterior expected gDNA counts
   from step 6 to update `A(x)`.

Stop after one refinement unless benchmarks show instability. The model should
be regularized enough that it does not chase expression or random coverage
spikes.

## Diagnostics To Emit

Add a calibration block to `summary.json`:

```text
capture_exposure:
  mode: uniform | learned | bed
  on_off_ratio: ...
  effective_genome_fraction: ...
  n_tiles: ...
  n_on_tiles: ...
  heldout_loglik_gain: ...
```

Add per-locus diagnostics to `loci.feather`:

- `gdna_eff_len_unweighted`
- `gdna_eff_len_weighted`
- `capture_weighted_fraction = weighted / unweighted`
- `capture_on_target_bp_proxy`
- `capture_off_target_bp_proxy`
- `capture_model_mode`

For annotated BAM debugging, add pool-level runner-up tags before changing hard
assignment policy:

- best RNA-vs-gDNA log posterior margin;
- gDNA posterior summed across the gDNA component;
- RNA posterior summed across mRNA+nRNA components.

These tags are essential for verifying that exposure weighting is fixing the
actual posterior margin, not only changing the winning component label.

## Acceptance Tests

### Synthetic Unit Test

Construct one 100 kb locus with a 10 kb target island and 90 kb off-target
sequence. Simulate gDNA with `A_on/A_off = 100`. Assert:

```text
weighted_gdna_eff_len ≈ 10 kb + 90 kb / 100
unweighted_gdna_eff_len ≈ 100 kb
```

Then simulate ambiguous unspliced fragments in the target island and verify that
the weighted model increases gDNA posterior relative to the unweighted model.

### Non-Capture Regression

Simulate uniform gDNA with `A(x)=1`. The learned model must classify as uniform,
and `gdna_eff_len_weighted / gdna_eff_len_unweighted` should be within a tight
tolerance, for example `[0.95, 1.05]` for normal loci.

### Real Hybrid-Capture Benchmark

Use the VCaP RNA/gDNA mix. Primary metric: gDNA-source fragments called RNA in
the no-multimap run, stratified by ambiguous unspliced category and by the five
hotspots in `docs/benchmarks/vcap_no_mm_deep_regions_2026-05-18.md`.

Expected direction:

- large drop in gDNA -> RNA calls in AR, ZNF462, MALAT1/TALAM1, FLG2, and the
  chr2 long-nRNA region;
- no large loss of true RNA recall for spliced mature-RNA fragments;
- `capture_weighted_fraction` much less than 1 in large capture-shaped loci;
- near-unity weighted fraction in non-capture synthetic benchmarks.

## Open Design Questions

1. Should the first implementation add `log W(f)` to gDNA scoring immediately,
   or start by weighting only `gdna_eff_len` and priors? The full model requires
   both; an effective-length-only prototype is easier but theoretically
   incomplete.
2. How aggressively should ambiguous exon-contained unspliced reads influence
   no-BED exposure learning after the first EM pass? A conservative cap is
   advisable.
3. Should mature RNA and nRNA effective lengths be capture-weighted in the same
   milestone or a follow-up? The gDNA fix is the urgent failure mode, but the
   full capture model should eventually weight all molecule classes.
4. What is the best tile scale? Start with 1 kb fixed windows intersected with
   `region_df`, then tune using held-out likelihood and hotspot correction.

## Recommended Path

Implement this in two stages.

Stage 1 should be an exposure-weighted gDNA model with an optional BED path and
a conservative no-BED learner. It should produce weighted gDNA effective length,
weighted gDNA priors, `log W(f)` scoring, and diagnostics, while preserving the
uniform fast path.

Stage 2 should generalize the same exposure field to mature RNA and nRNA
effective lengths and add more robust iterative posterior-based exposure
learning. This is where Rigel becomes capture-aware as a quantifier, not only as
a gDNA/RNA classifier.