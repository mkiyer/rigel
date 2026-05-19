# gDNA Exposure Model v2

Status: implementation plan draft, not implemented

Date: 2026-05-18

## Goal

Fix hybrid-capture gDNA under-calling with the smallest model that is still
theoretically coherent.

The v1 plan introduced a general genomic exposure field `A(x)`. That is the
right umbrella, but too broad for a first patch. The v2 implementation should
start with a two-state capture model:

```text
on-target weight  = 1
off-target weight = lambda = 1 / c
```

where `c` is a sample-level capture enrichment ratio learned during
calibration. For uniform polyA/ribo-minus libraries, `c ~= 1` and Rigel uses the
current unweighted path exactly. For capture libraries, `c >> 1` and the gDNA
opportunity space shrinks toward the on-target footprint.

nRNA remains unspliced by definition. Spliced alignments are mature RNA evidence
or alignment artifacts, never nRNA evidence.

## What We Are Taking From The Alternative Proposal

The other proposal has several good simplifications worth adopting.

1. The equivalent-uniform-length identity is the right first implementation:

   ```text
   L_gdna_capture = L_on + L_off / c
   ```

   This is just a two-state exposure field normalized to on-target density.

2. A single sample-level `c` is a good starting point. Per-locus enrichment
   ratios are identifiable only with strong local evidence and should be
   deferred.

3. The implementation should be staged. First prove that capture-aware gDNA
   effective length fixes the VCaP false-RNA hotspots; then add the per-fragment
   positional term before considering the model production-ready.

4. Optional BED support should share the same internal model. A BED file should
   provide target geometry, not a separate algorithm.

5. The model should report a library diagnostic rather than requiring users to
   choose a capture mode.

## What We Should Not Adopt As-Is

1. Non-expressed exons are not guaranteed pure gDNA evidence. They are a useful
   proxy in exon-capture assays, but they can still contain unannotated RNA,
   nRNA, alignment artifacts, or panel-specific off-target behavior. Use them
   with robust filters, not as a single source of truth.

2. A per-locus HMM is too much for the first implementation. Start with target
   intervals from BED or a simple target proxy, plus global `c`.

3. Per-locus `c` should be deferred. It is tempting, but it risks turning the
   exposure model into an arbitrary local density smoother.

4. EM posterior-driven exposure refinement should be deferred until the
   denominator-only hypothesis has been benchmarked. It introduces circularity
   unless cross-fit carefully.

## Model

For a `MultiLocus` `M`, split the genomic opportunity into target and
non-target start positions. With on-target weight 1 and off-target weight
`lambda = 1 / c`:

```text
P(f | gDNA_M) ∝ h_G(ell_f) q_G(f) w(f) / Ltilde_gDNA(M)

Ltilde_gDNA(M) = sum_ell h_G(ell) [N_on(M, ell) + lambda N_off(M, ell)]
```

`N_on + lambda N_off` is shorthand for summing the exposure weight of each
valid fragment start. In the simplest implementation, a start is classified by
the fragment midpoint. A more exact implementation integrates the two-state
target weight across the fragment footprint. Either way, this is the weighted
version of `_exposure.gdna_eff_len_for_loci()`.

For a one-interval locus far from contig boundaries, the intuition is:

```text
Ltilde_gDNA ~= L_on + lambda L_off + E[ell] - 1 boundary terms
```

The boundary terms should be handled by the same FL-marginal start-window logic
already used for gDNA overlap effective length, not by a hand-written scalar
span formula.

The per-fragment term is:

```text
w(f) = 1       if the fragment footprint is on-target
w(f) = lambda  if the fragment footprint is off-target
```

For fragments partially overlapping target intervals, use the same convention
as the effective-length projector: midpoint in the first implementation, or
footprint-average exposure in the exact implementation.

Adding `log w(f)` to the gDNA score and using weighted `Ltilde_gDNA` makes the
likelihood mass-conserving. A denominator-only prototype is acceptable for
debugging, but production should use both numerator and denominator weights.

## Target Geometry

Represent target geometry as a compact interval model:

```text
GdnaExposureModel:
  mode: uniform | bed | learned_annotation | learned_coverage
  c: float
  lambda: float
  target_intervals_by_ref: sorted half-open intervals
```

If `mode == uniform`, target intervals are ignored and `lambda = 1`.

### Source Priority

1. **Explicit BED**: if provided, target intervals come from the BED, optionally
   expanded by a small shoulder tied to the gDNA fragment-length distribution.
   Calibration still estimates `c` from the BAM.
2. **Learned simple target map**: no BED. Use a conservative automatic map from
   annotation plus coverage proxies.
3. **Uniform fallback**: if the capture diagnostic is weak or evidence is
   insufficient, set `c = 1` and use the current model.

### No-BED Target Map v2

The first no-BED implementation should be deliberately simple:

1. Build a regional exposure table from the existing calibration region
   partition. This uses the machinery Rigel already has: intergenic contained
   counts, intron contained counts, exon-intron boundary counts, FL-aware
   exposure, and NB shrinkage toward global density channels.
2. Use annotated exonic intervals and eligible exon-intron shoulders as a weak
   target prior when regional evidence is sparse. This matches common exon/panel
   capture behavior and is cheap because the index already has exon/intron
   partitions.
3. Optionally refine with high-density gDNA-proxy tiles after the global `c`
   estimate is available. Use fixed 1 kb tiles intersected with `region_df`.
4. Do not use mature spliced RNA coverage as positive target evidence.
5. Do not use EM posterior labels in v2 target learning.

This will not perfectly recover every panel, but it is straightforward and
should address the observed failure mode: large transcript/nRNA spans with a
small capture-accessible footprint.

## Regional-Density Exposure Weights

The existing calibration region partition is a good substrate for learning
capture exposure. The correct object to learn is not a normalized mass fraction;
it is a relative sampling density.

For each region `r`, define:

```text
E_r       = FL-aware opportunity/exposure for the region
Y_r       = conservative gDNA-proxy count for the region
rho_r     = shrunk regional gDNA density estimate
A_r       = relative exposure weight, proportional to rho_r
```

The current code already has most ingredients during calibration:

- `CalibrationScanPayload.per_region_counts` for intergenic and intronic
  contained counts;
- `u_left` / `u_right` for exon-intron boundary counts;
- global density channels and `kappa` overdispersion estimates;
- `shrink_to_loco(...)`, which gives the right empirical-Bayes shape for
  shrinking small-exposure estimates toward a global density.

### Correct Use

Use regional density to construct an exposure field:

```text
A_r = clamp(shrink(rho_r) / rho_ref, lambda_floor, 1)
weighted_exposure = sum_r A_r * E_r
```

`rho_ref` can be a high robust quantile of `rho_r` in capture libraries, or the
global density in uniform libraries. The absolute scale is arbitrary as long as
the same `A_r` is used in both gDNA score and gDNA effective length.

This regional table can feed either implementation style:

- two-state: threshold `A_r` into target/off-target and estimate `c`;
- continuous: keep `A_r` directly as a piecewise-constant exposure field.

The two-state path is easier to reason about. The continuous path may be more
faithful and still fits the same weighted effective-length code.

### Incorrect Use

Do not assign each region a weight equal to its fraction of total global gDNA
mass and then multiply the region's length by that weight.

If `m_r = rho_r * E_r` and `W_r = m_r / sum(m)`, then a weighted denominator
like `sum_r W_r * E_r` becomes:

```text
sum_r rho_r * E_r^2 / sum(m)
```

That weights by length twice and is not the gDNA opportunity space. It also
changes when the region partition is split into smaller intervals, which a
physical likelihood must not do.

The invariant is:

```text
expected mass in region r = rho_ref * A_r * E_r
```

Therefore `A_r` should be proportional to density, not to total mass fraction.

### Uncertainty

Uncertainty should enter through shrinkage, not through ad hoc powers like
`rho_r ** gamma` as the main model.

Small regions and exon-intron boundaries have high variance. Their regional
densities should shrink strongly toward the appropriate global density channel
unless enough local evidence supports a high-exposure state. Larger regions with
many proxy fragments can keep more local signal.

A later implementation can carry posterior variance or credible intervals for
`A_r`, but v2 should start with one shrunk posterior mean per region plus
diagnostics:

- `region_gdna_density_raw`
- `region_gdna_density_shrunk`
- `region_exposure_weight`
- `region_exposure_uncertainty_proxy`

This is enough to avoid the worst instability while keeping the model simple.

## Estimating The Enrichment Ratio `c`

Estimate one robust sample-level capture ratio during calibration.

### With BED

Compute density in target and non-target proxy regions:

```text
c = rho_on / rho_off
```

Use gDNA-proxy counts only: intergenic-contained, intron-contained,
exon-intron boundary events, opposite-strand unspliced in stranded libraries,
and exonic unspliced counts from transcripts with no mature spliced support.
The exonic non-expressed proxy is useful here, but it should be winsorized and
filtered rather than trusted blindly.

### Without BED

Use a robust density-spread estimator over proxy tiles:

```text
d_i = proxy_gDNA_count_i / FL_exposure_i
capture_index = Q95(log d_i) - Q25(log d_i)
c = exp(capture_index)
```

Apply minimum support checks and caps:

```text
if n_proxy_tiles < threshold: c = 1
if c < c_min_capture: c = 1
c = min(c, c_max)
```

This is simpler than a mixture model and likely enough to separate uniform
libraries from strong capture libraries. A two-component mixture can be added
later if quantile estimates are unstable.

Recommended initial constants for experiments:

```text
c_min_capture = 3
c_max = 200
tile_size = 1000 bp
lambda_floor = 1 / c_max
```

These should be config values but not exposed as prominent CLI knobs initially.

## Implementation Stages

### Stage 0: Diagnostics Only

Add a debug script or calibration diagnostic that estimates `c`, target
fraction, and weighted/unweighted gDNA effective length for existing outputs.
Run it on the VCaP no-multimap BAM and the five hotspot regions.

This script should first build the regional-density exposure table described
above and compare three candidate denominators:

```text
unweighted span / current gdna_eff_len
two-state target/off-target weighted length
continuous regional-density weighted length
```

This stage tests the hypothesis before touching the EM path.

### Stage 1: Denominator And Prior Prototype

Add Python-side weighted gDNA exposure for EM component metadata:

- `weighted_gdna_eff_len_for_loci(...)`
- `RegionGdnaExposureTable` or equivalent compact interval-weight model
- weighted `expected_gdna_count_global(...)`
- `loci.feather` diagnostics:
  - `gdna_eff_len_unweighted`
  - `gdna_eff_len_weighted`
  - `gdna_eff_len_capture_ratio`
  - `capture_c`
  - `capture_target_bp_proxy`
  - `capture_offtarget_bp_proxy`

This stage plugs into the existing `gdna_eff_len` and `gdna_prior_count` arrays
with minimal native changes. It is theoretically incomplete because it does not
yet add `log w(f)` to each gDNA fragment score, but it is the fastest way to
measure whether denominator inflation is the dominant VCaP failure.

### Stage 2: Per-Fragment gDNA Weight

Add the numerator term:

```text
gdna_log_lik[f] += log w(f)
```

The cleanest place is native scoring, where fragment coordinates are already
available. If passing interval state into C++ is too invasive, compute a
per-unit `gdna_exposure_log_weight` in the Python routing layer as an interim
step, but the final model should score it in the same pass as gDNA likelihoods.

After Stage 2, the model is mass-conserving and ready for production testing.

### Stage 3: Optional Coverage Refinement

Only if Stage 1/2 under-correct or over-correct, refine the no-BED target map
with raw pre-EM coverage segmentation:

- fixed 1 kb tiles;
- conservative gDNA-proxy counts;
- shrink high-density calls toward the annotation target prior at low support;
- no EM posterior feedback in the first refinement.

Defer HMMs, per-locus `c`, and iterative posterior-based learning until there
is evidence that the simple model is insufficient.

## Safeguards

1. **Uniform exact fallback**: when `c = 1`, weighted effective lengths must be
   exactly or numerically identical to current unweighted values.
2. **No expression-as-target shortcut**: spliced mature RNA coverage cannot be
   used to mark targets.
3. **No circular learning in v2**: EM posterior assignments cannot define the
   target mask used by that same EM run.
4. **Low-complexity exposure**: one global `c`; target/off-target intervals;
   no free per-locus density unless added in a later empirical-Bayes model.
5. **nRNA invariant**: nRNA candidates should receive only unspliced evidence.
   Any spliced+nRNA annotations in debug output should be audited separately.

## Expected Effect On VCaP Hotspots

The five no-multimap hotspots are exactly the pattern this model targets:

- large `MultiLocus` or long nRNA span;
- small capture-accessible footprint;
- many `NH=1`, low-`NM`, `ZS=unspliced` gDNA-source fragments;
- gDNA denominator inflated by off-target span.

For a 100 kb locus with 10 kb target and `c = 100`:

```text
unweighted denominator ~= 100 kb
weighted denominator   ~= 10 kb + 90 kb / 100 = 10.9 kb
gDNA log posterior gain ~= log(100 / 10.9) = 2.2 nats
```

For a 1.6 Mb locus with 8 kb target and `c = 100`:

```text
weighted denominator ~= 8 kb + 1592 kb / 100 = 23.9 kb
gain ~= log(1600 / 23.9) = 4.2 nats
```

That is large enough to change ambiguous unspliced posterior competition
without altering the mature spliced RNA model.

## Tests

1. **Geometry unit test**: one locus, target/off-target split, delta fragment
   length. Assert weighted length equals `L_on + L_off / c` plus the expected
   FL-overlap boundary convention.
2. **Uniform regression**: `c = 1` gives the same `gdna_eff_len` and EM results
   as the current code.
3. **Synthetic capture test**: simulate on/off gDNA coverage with known `c` and
   verify the estimator recovers the ratio within a broad tolerance.
4. **Regional-weight invariance test**: splitting one region into two adjacent
   regions with the same density must not change weighted gDNA effective length.
   This guards against accidental mass-fraction weighting.
5. **VCaP hotspot benchmark**: rerun no-multimap VCaP and report gDNA-source
   fragments called RNA in the five hotspot windows from
   `docs/benchmarks/vcap_no_mm_deep_regions_2026-05-18.md`.
6. **Non-capture benchmark**: rerun existing synthetic/polyA-style benchmarks
   and require near-zero change when the capture diagnostic chooses `c = 1`.

## Open Questions

1. Is annotation exonic target prior sufficient for this VCaP capture assay, or
   do we need the Stage 3 coverage refinement immediately?
2. Should Stage 1 denominator-only be an internal experiment only, or should it
   be hidden behind an explicit experimental config flag?
3. How should target shoulders be defined: fixed bp, mean gDNA FL, or a small
   quantile of the gDNA FL distribution?
4. After gDNA is fixed, do mature RNA and nRNA effective lengths need the same
   target/off-target exposure correction for unbiased quantification? This is
   likely yes, but it should be a follow-up after the gDNA false-RNA mode is
   controlled.

## Recommended Next Step

Implement Stage 0 as a debug/profiling script first. It should build regional
gDNA density/exposure weights from the existing calibration payload, then compute
`c`, `L_on`, `L_off`, unweighted `gdna_eff_len`, two-state weighted
`gdna_eff_len`, and continuous regional-density weighted `gdna_eff_len` for the
five VCaP hotspots. If the predicted gDNA log-posterior gains match the observed
failure scale, proceed to Stage 1 in the calibration code.