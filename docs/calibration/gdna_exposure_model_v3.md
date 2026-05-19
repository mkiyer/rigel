# gDNA Exposure Model v3

Status: implementation plan draft, not implemented

Date: 2026-05-18

## Decision

Pursue a region-partition exposure model.

The cleanest first implementation is not a separate target detector and not a
free-form local smoother. It is a compact relative exposure field learned from
Rigel's existing calibration region partition.

For each genomic region `r`, estimate a shrunk gDNA proxy density and convert it
to a relative exposure weight:

```text
E_r       = FL-aware opportunity/exposure for region r
Y_r       = conservative gDNA-proxy count in region r
rho_r     = EB-shrunk regional gDNA density
A_r       = clamp(rho_r / rho_ref, lambda_floor, 1)
```

Then use `A_r` in both places where gDNA exposure matters:

```text
gDNA score:          log P(f | gDNA) += log A(fragment)
gDNA effective len:  L_g = sum valid starts A(start, ell)
```

This is theoretically the same as the earlier capture exposure field, but it is
implemented with infrastructure Rigel already has.

## Why This Is The Right Direction

Hybrid capture creates spatially uneven gDNA sampling. A large `MultiLocus` can
span mostly depleted intron/off-target sequence while the fragments pile up on a
small captured footprint. Normalizing gDNA by the whole unweighted footprint
penalizes gDNA and routes ambiguous unspliced fragments to mRNA or nRNA.

The region partition is already the right grid for measuring this. It separates:

- intergenic contained gDNA evidence;
- intronic contained gDNA evidence;
- exon-intron boundary gDNA evidence;
- exon-contained gDNA projected from boundary evidence.

For polyA/ribo-minus libraries, the shrunk regional densities should be nearly
flat, so `A_r ~= 1` and the model becomes the current uniform model. For capture
libraries, high-density regions keep `A_r` near 1 and depleted regions get small
weights.

## Important Correction

Do not weight a region by its fraction of total gDNA mass.

If `m_r = rho_r * E_r` and `W_r = m_r / sum(m)`, then a denominator like
`sum_r W_r * E_r` becomes:

```text
sum_r rho_r * E_r^2 / sum(m)
```

That weights by region length twice and changes when the same physical interval
is split into smaller regions. A physical likelihood must be invariant to such a
partition split.

The correct quantity is relative density:

```text
weighted exposure = sum_r A_r * E_r
A_r proportional to rho_r, not rho_r * E_r
```

This keeps the invariant:

```text
expected gDNA mass in region r = rho_ref * A_r * E_r
```

## Minimal Model

### Region Weights

Use one piecewise-constant weight per calibration region:

```text
RegionGdnaExposureTable:
  region_id
  ref_id, start, end
  region_type
  exposure_bp
  proxy_count
  rho_raw
  rho_shrunk
  exposure_weight
  uncertainty_proxy
```

`exposure_weight` is dimensionless. It is normalized so that a typical captured
or high-density region has weight near 1.

Recommended reference density:

```text
rho_ref = robust high quantile of shrunk regional densities, e.g. Q95
```

Clamp weights:

```text
A_r = clamp(rho_shrunk_r / rho_ref, lambda_floor, 1)
```

Initial values for experiments:

```text
lambda_floor = 0.005       # equivalent to max enrichment c = 200
min_regions_for_model = configurable diagnostic threshold
```

If the regional density distribution is not meaningfully overdispersed, set all
weights to 1 and mark the exposure model as `uniform`.

### Conservative gDNA Proxy Counts

Only use fragment categories that are gDNA-dominated before EM feedback:

- intergenic-contained counts;
- intron-contained counts, strand-corrected when the strand model is strong;
- exon-intron boundary counts, strand-corrected when possible;
- exon-contained unspliced counts only in conservative situations, such as
  opposite-strand evidence in stranded libraries or transcripts with no mature
  spliced support.

For v3, the safest implementation is to start with the categories already used
by calibration and avoid adding exon-contained ambiguous counts until the first
diagnostic run shows they are necessary.

### Exon Regions

Exon-contained fragments are exactly where gDNA and mature RNA are hardest to
separate. Do not estimate exon gDNA exposure directly from total exon coverage.

Instead, assign exon-region gDNA density from boundary evidence:

```text
rho_exon_region = shrunk density from eligible exon-intron boundary sides
```

If a small exon has no eligible boundary evidence, shrink to the exon-intron
global density or to the local neighboring boundary density. This is not perfect,
but it prevents mature RNA expression from becoming target evidence.

## Uncertainty

Uncertainty should be modeled through shrinkage first, not through powers of
density and not through an immediate new EM prior-strength system.

Use the existing empirical-Bayes shape:

```text
rho_shrunk = (Y_r + kappa_class * rho_global_class) / (E_r + kappa_class)
```

This is equivalent to a Gamma-Poisson style posterior mean and has the behavior
we want:

- tiny regions shrink strongly to the class/global density;
- large high-count regions retain local signal;
- exon-intron boundary regions do not create noisy spikes from one or two
  fragments.

Emit an uncertainty proxy for diagnostics, for example posterior coefficient of
variation or effective precision:

```text
precision_proxy = E_r + kappa_class
uncertainty_proxy = 1 / sqrt(max(precision_proxy, eps))
```

Do not use this uncertainty to change EM prior strength in the first production
patch. Weighted expected gDNA counts already alter the gDNA prior. Propagating
variance into Dirichlet strength is a good follow-up, but it is not needed to
test or fix the main denominator failure.

## Weighted Effective Length

Add a weighted sibling of the current gDNA overlap effective length:

```text
weighted_gdna_eff_len_for_loci(loci, ref_lengths, fl, exposure_table)
```

Use one convention for both denominator and per-fragment score. The v3 default
should be the midpoint convention because it is simple and fast:

```text
for a start s and length ell:
  midpoint = s + floor(ell / 2)
  A(start, ell) = exposure_weight_at(midpoint)
```

Then:

```text
L_g(M) = sum_ell h_G(ell) sum_valid_starts A(start, ell)
```

This remains compatible with the existing start-window logic in
`gdna_eff_len_for_loci()`. The exact footprint-average version can be added
later if midpoint classification fails near sharp target boundaries.

Required invariant:

```text
If all A_r = 1, weighted_gdna_eff_len_for_loci == gdna_eff_len_for_loci.
```

## Per-Fragment Score

The same regional weight must enter the gDNA numerator:

```text
gdna_log_lik[f] += log(max(A(fragment_midpoint), lambda_floor))
```

The denominator-only version is useful as a diagnostic experiment, but it should
not be treated as the final model. Production should include both weighted
effective length and per-fragment gDNA weight.

## Priors

Update `expected_gdna_count_global()` to project global densities onto weighted
region exposure:

```text
eta_g = sum_r rho_global_class(r) * A_r * E_r
```

For exon regions, use the boundary-derived exon density channel and weighted
exon-contained exposure. For boundary-crossing expected mass, weight the
boundary side by the adjacent exon region's exposure weight.

Keep the current Bayesian prior structure otherwise. Do not introduce a new
variance-weighted Dirichlet strength until we have benchmark evidence that the
simple weighted expected-count prior is insufficient.

## Optional BED

A target BED should be optional and should not create a separate algorithm.

If provided, it can be used in two conservative ways:

1. as a prior or diagnostic label for regional exposure states;
2. as a fallback target geometry when regional evidence is sparse.

Even with BED, estimate the enrichment/weights from the BAM. BED geometry is not
always assay-exact, and data-driven weights protect against stale target files.

## Implementation Stages

### Stage 0: Diagnostic Script

Add a script under `scripts/debug/` that builds a regional exposure table and
evaluates the five VCaP hotspots.

For each hotspot report:

- current unweighted `gdna_eff_len`;
- weighted regional-density `gdna_eff_len`;
- predicted gDNA log-posterior gain `log(unweighted / weighted)`;
- top contributing high-weight and low-weight regions;
- uncertainty summaries for the regions contributing most exposure.

This can be run before changing production code. If the predicted gains are not
large enough in the known false-RNA regions, stop and investigate before
touching the EM path.

### Stage 1: Python Exposure Table And Weighted Denominator

Add a production `RegionGdnaExposureTable` builder during calibration. Thread it
into prior assembly and compute weighted gDNA effective lengths.

Outputs:

- `summary.json` exposure diagnostics:
  - `mode`: `uniform` or `regional`
  - `capture_index`
  - `rho_ref`
  - `lambda_floor`
  - number of regions with weight below thresholds
- `loci.feather` diagnostics:
  - `gdna_eff_len_unweighted`
  - `gdna_eff_len_weighted`
  - `gdna_eff_len_weight_ratio`
  - regional exposure summary statistics

Run this behind an experimental config flag until VCaP and non-capture
benchmarks agree with expectations.

### Stage 2: Matching gDNA Score Weight

Add the per-fragment numerator term in the scorer or routing layer. This makes
the model mass-conserving.

Preferred final location: native scoring, because fragment coordinates and gDNA
likelihood are already there. An interim Python-side per-unit log-weight is
acceptable for experiments if it avoids a larger C++ interface change.

### Stage 3: Only If Needed

Defer these until Stage 0-2 results show a real need:

- footprint-average weights instead of midpoint weights;
- per-locus enrichment ratios;
- HMM or change-point target segmentation;
- EM-posterior-driven exposure refinement;
- using exposure uncertainty to modulate Dirichlet prior strength;
- applying the same exposure weights to mature RNA and nRNA effective lengths.

The last item is likely necessary for fully unbiased capture quantification, but
it should not block the first gDNA false-RNA fix.

## Tests

1. **Uniform identity**: all weights 1 reproduces current gDNA effective length,
   priors, and EM output within numerical tolerance.
2. **Partition invariance**: splitting a region into adjacent subregions with
   identical density and weight does not change weighted effective length.
3. **Synthetic regional capture**: simulate a locus with high-density exon
   regions and low-density introns. Verify the weighted denominator matches the
   known exposure and increases gDNA posterior for ambiguous unspliced gDNA
   fragments.
4. **Uncertainty shrinkage**: tiny regions with one fragment shrink strongly;
   large regions with matching density retain signal.
5. **VCaP hotspot benchmark**: rerun the no-multimap VCaP case and measure
   gDNA-source fragments called RNA in the five documented hotspots.
6. **Non-capture benchmark**: existing synthetic/polyA-style benchmarks should
   choose `uniform` or near-uniform weights and show near-zero quantification
   drift.

## Recommendation

Implement the regional exposure diagnostic first. If it predicts meaningful
gDNA log-posterior gains in the known VCaP false-RNA hotspots, proceed with the
simple regional model in production.

Keep v3 small:

- regional density weights;
- EB shrinkage;
- weighted gDNA effective length;
- matching per-fragment gDNA weight;
- diagnostics.

Do not add per-locus `c`, HMMs, posterior feedback, covariance propagation, or
RNA/nRNA exposure weighting until this simple model has been benchmarked.