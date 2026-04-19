# Unified Joint-Likelihood gDNA Calibration

_Date: 2026-04-18 · Author: Copilot analysis, reviewed by MK_

## Motivation

The current calibration architecture has two independent pathways
(strand-based and density-based) blended by a heuristic weight.
End-to-end validation on unstranded + high-gDNA STAR-aligned data
([stress_test_findings.md](stress_test_findings.md)) exposed two
structural problems:

1. **The P₁₀ density estimator collapses to zero** under realistic
   sparsity. Even at 50% gDNA contamination the per-bp fragment rate
   is only λ ≈ 10⁻³, so the majority of 500 bp regions are expected
   to have zero reads under the null (see §1 below). P₁₀ is not a
   robust estimator in this regime.

2. **Unmapped reads are ignored.** STAR flags reads that multimap to
   too many loci as unmapped (`uT:3`). These are overwhelmingly gDNA
   from repetitive regions and represent 3–15% of a typical library —
   an enormous information source for λ_G that we currently discard.

This document proposes replacing the two-pathway blend with a single
joint log-likelihood over three information sources, eliminating all
arbitrary cutoffs (P₁₀, w_strand blend weight) and rescuing
unmapped-read signal.

---

## 1. The sparsity problem (why P₁₀ must go)

### 1.1. Expected zeros under the null

At $N_\text{gDNA}$ fragments over a genome of size $G$, the per-bp
fragment-start rate is

$$\lambda_G = \frac{N_\text{gDNA}}{G}$$

For a region of genomic span $S$ bp, the expected count is
$\mu = \lambda_G \cdot S$ and the probability of observing zero reads
is

$$P(k = 0) = e^{-\lambda_G S}$$

Concrete example: 5 M gDNA fragments, $G = 4 \times 10^9$ bp,
$\lambda_G = 1.25 \times 10^{-3}$, region span $S = 500$ bp:

$$\mu = 0.625, \quad P(k=0) = e^{-0.625} \approx 0.535$$

More than half of all pure-gDNA regions are expected to have zero
reads. For a P₁₀ to be nonzero, we need $P(k=0) < 0.10$, which
requires $\mu > \ln 10 \approx 2.3$ — i.e. regions of ~1,800 bp or
3–4× more sequencing depth. Neither is realistic for a general-purpose
tool.

### 1.2. Scaling with depth

| Total reads | gDNA fraction | $\lambda_G$ | $\mu$ (500 bp) | $P(k=0)$ | P₁₀ viable? |
|---|---|---|---|---|---|
| 10 M | 50% | 1.25e-3 | 0.625 | 53.5% | No |
| 20 M | 50% | 2.50e-3 | 1.250 | 28.7% | Marginal |
| 50 M | 50% | 6.25e-3 | 3.125 | 4.4% | Yes |
| 100 M | 50% | 1.25e-2 | 6.250 | 0.2% | Yes |
| 20 M | 10% | 5.00e-4 | 0.250 | 77.9% | No |

The P₁₀ estimator is fundamentally incompatible with the
low-to-moderate depth, moderate-gDNA regime that characterises most
public RNA-seq datasets (10–30 M reads, gDNA fraction 5–30%).

### 1.3. What P₁₀ was trying to do

The 10th percentile of per-region density was a proxy for "the
RNA-free tail of the density distribution." The underlying idea is
sound: gDNA is approximately uniform; RNA is concentrated in exons;
the lowest-density regions should reflect gDNA alone. But a hard
percentile cutoff is the wrong statistic when the null distribution
has most of its mass at zero.

---

## 2. Candidate replacements for P₁₀

### 2.1. ★★★ Poisson mixture MLE (recommended primary estimator)

Model each region $i$ as a two-component Poisson mixture:

$$k_i \sim (1 - \pi_i)\,\text{Poisson}(\lambda_G\, f_i\, S_i) \;+\; \pi_i\,\text{Poisson}(\lambda_G\, f_i\, S_i + \mu_{R,i}\, S_i)$$

where:

- $k_i$ = observed unspliced fragment count in region $i$
- $S_i$ = genomic span (bp)
- $f_i \in [0,1]$ = mappable fraction of region $i$ (continuous weight
  derived from the mappable BED, replacing the hard in/out mask)
- $\lambda_G$ = per-bp gDNA fragment rate (the quantity we want)
- $\pi_i$ = posterior probability that region $i$ contains RNA signal
- $\mu_{R,i}$ = per-bp RNA rate in region $i$ (nuisance)

The EM updates are standard two-component Poisson mixture:

1. **E-step:** compute $\pi_i$ for each region given current
   $\hat\lambda_G$ and $\hat\mu_R$
2. **M-step:** update $\hat\lambda_G$ as the
   posterior-weighted-average rate across regions assigned to the null
   component, weighted by $f_i S_i$

Properties:

- **Zeros are informative.** Under the null component, a region with
  $k_i = 0$ contributes $e^{-\lambda_G f_i S_i}$ to the likelihood —
  its span acts as exposure. Many zeros with large $f_i S_i$ directly
  constrain λ_G from above.
- **No arbitrary percentile cutoff.** The mixture handles the boundary
  between "gDNA-only" and "gDNA+RNA" regions probabilistically.
- **Partial mappability as a weight.** The $f_i$ term eliminates the
  hard mask decision. Regions with $f_i = 0.3$ contribute 30% of their
  span to the effective exposure — no region is thrown away.
- **Uncertainty is free.** Observed Fisher information at
  $\hat\lambda_G$ gives a confidence interval for calibration
  diagnostics.

Convergence is fast because λ_G is 1-dimensional and the profile
likelihood is log-concave.

### 2.2. ★★★ Zero-count moment estimator (recommended bound-check)

Under the null, the expected fraction of zero-count regions is

$$\mathbb{E}[\text{frac zeros}] = \frac{1}{N}\sum_i e^{-\lambda_G f_i S_i}$$

RNA contamination can only push a region *out* of the zero class (by
adding counts), so the observed zero fraction is an **upper bound** on
the null zero fraction. This gives a **lower bound** on λ_G:

$$\hat\lambda_G^{\text{lb}} = \arg\min_\lambda \left|\sum_i \mathbf{1}[k_i = 0] - \sum_i e^{-\lambda f_i S_i}\right|$$

This is a 1-D root-finding problem (bisection) and uses only the
zero-indicator — completely robust to the RNA tail. Report as a
diagnostic / sanity check alongside the mixture MLE.

### 2.3. ★★ Trimmed-sum with data-adaptive cutoff

Rank regions by $k_i / (f_i S_i)$. Iteratively fit a Poisson MLE on
the lower partition, then expand/shrink the partition until the
included regions pass a Poisson goodness-of-fit test (KS or
chi-squared). This is the data-adaptive generalisation of P₁₀ — the
"10" is replaced by whatever fraction fits the null.

Disadvantage: still discards part of the data; less principled than the
mixture model.

### 2.4. ★ Intergenic-only pathway

Restrict the MLE to regions the rigel index labels as intergenic.
Simple and fast, but vulnerable to unannotated transcription in
intergenic space.

### 2.5. Recommendation

Use **§2.1 (Poisson mixture MLE)** as the primary density-pathway
estimator. Use **§2.2 (zero-count moment estimator)** as a diagnostic
lower bound reported in `summary.json`. Retire P₁₀ entirely.

---

## 3. Rescuing unmapped reads

### 3.1. The information source

STAR writes a `uT` BAM tag on unmapped reads indicating the reason for
non-mapping:

| `uT` value | STAR reason | Dominant biological source | Signal value |
|---|---|---|---|
| 0 | No acceptable seed | Foreign organism / adapter dimer | None (discard) |
| 1 | Too short | Adapter / artefact | None (discard) |
| 2 | Too many mismatches | Divergent sequence, fusion breakpoints, foreign DNA | Weak / ambiguous |
| 3 | Too many multimapping loci | Repetitive-element gDNA (LINE/SINE/SD) + repeat-RNA | **High (primary gDNA signal)** |
| 4 | Unmapped mate of mapped read | Mixed | Low |

The `uT:3` category — reads that exceeded `--outFilterMultimapNmax` —
is the goldmine. These reads map to so many loci that STAR refuses to
report them. In a gDNA-contaminated library, most of these come from
high-copy-number repeats that are uniformly sampled by gDNA
fragmentation but too repetitive for unique alignment.

### 3.2. The denominator problem

To convert a count $N_\text{uT3}$ into a density, we need a genomic
span. The natural denominator is the **total span of the genome that is
unmappable at a given read length** — the inverse of the mappable BED.

Available resource: `alignable_unmappable_rl125_grch38_star.bed`
provides the truly unmappable regions for 125 bp reads under STAR
defaults. Let $S_\text{rep}$ denote the total span of this BED file.

### 3.3. Confounders in the uT:3 compartment

Not all uT:3 reads are gDNA. Two biological sources of repeat-derived
RNA exist:

1. **Active transposable elements (LINE-1, Alu, SVA).** These are
   transcribed and can be abundant (0.5–3% of total RNA in most
   tissues; higher in testes, placenta, some cancers). They generate
   reads that are genuinely unmappable.

2. **Satellite / pericentromeric transcription.** Low-level but real.

These are real RNA molecules landing in the uT:3 bin. We model this
with a prior RNA fraction:

$$\pi_\text{RNA} \approx 0.05 \quad \text{(default; Beta(2, 38) prior)}$$

meaning ~5% of uT:3 reads are assumed to be RNA. This prior can be
overridden per tissue type or refined from stranded libraries (where
repeat-RNA can be distinguished by strand asymmetry).

### 3.4. ★★★ Unmapped Poisson pathway (recommended)

Model:

$$N_\text{uT3,gDNA} = (1 - \pi_\text{RNA}) \cdot N_\text{uT3}$$

$$N_\text{uT3,gDNA} \sim \text{Poisson}(\lambda_G \cdot S_\text{rep})$$

Log-likelihood contribution:

$$\ell_\text{unmapped}(\lambda_G) = N_\text{uT3,gDNA}\,\log(\lambda_G\, S_\text{rep}) - \lambda_G\, S_\text{rep} - \log(N_\text{uT3,gDNA}!)$$

Key property: $S_\text{rep}$ is typically **200–400 Mb** — orders of
magnitude larger than any single region's span. This single term
carries enormous Fisher information:

$$I_\text{unmapped}(\lambda_G) = \frac{S_\text{rep}}{\lambda_G}$$

At $\lambda_G = 10^{-3}$ and $S_\text{rep} = 300\,\text{Mb}$:
$I = 3 \times 10^{11}$, giving a standard error of
$\sigma \approx 1.8 \times 10^{-6}$ — sub-percent precision from this
pathway alone.

Compare to the Fisher information from a single 500 bp mappable region:
$I = 500 / 10^{-3} = 5 \times 10^5$. The unmapped pathway has
**600,000×** more information per "observation" because the effective
exposure is the entire unmappable genome.

### 3.5. Sensitivity to π_RNA

Even at $\pi_\text{RNA} = 0.5$ (an extreme scenario), the point
estimate shifts by a factor of 2. But the likelihood is still sharply
peaked because $S_\text{rep}$ is so large. The uncertainty from
$\pi_\text{RNA}$ becomes the **dominant** source of calibration
uncertainty rather than sparsity — which is exactly where we want the
uncertainty to live, because $\pi_\text{RNA}$ is a biologically
meaningful quantity we can bound.

### 3.6. Practical implementation notes

- **BAM parsing:** The C++ scanner already reads unmapped records.
  Adding a `uT` tag counter is a one-line addition to `bam_scanner.cpp`
  (read the `uT` aux tag, increment a counter per value 0–4).
- **$S_\text{rep}$:** Computed once from the unmappable BED at index
  time and stored as a scalar in the index metadata. Read-length
  dependence is modest (125 bp vs 150 bp changes $S_\text{rep}$ by
  <5%); a single pre-computed value is sufficient.
- **$\pi_\text{RNA}$ prior:** Stored in `PipelineConfig` with a
  sensible default (0.05). Overridable via CLI
  (`--unmapped-rna-fraction`).

---

## 4. The strand pathway (retained, not replaced)

The existing strand-based calibration remains valid and powerful when
$\sigma \gg 0.5$. It provides a completely orthogonal signal: gDNA
generates opposite-strand fragments in exonic regions, which RNA does
not. No changes to its internal logic are needed.

Its contribution in the unified framework is a standard Poisson
log-likelihood on opposite-strand exonic counts:

$$\ell_\text{strand}(\lambda_G) = \log\text{Poisson}\!\left(N_\text{opp}\;\Big|\;\lambda_G \cdot S_\text{exon} \cdot \frac{1 - \sigma}{2}\right)$$

where $N_\text{opp}}$ is opposite-strand exonic fragment count,
$S_\text{exon}$ is total exonic span, and $\sigma$ is the estimated
strand specificity.

At $\sigma = 0.5$ (unstranded), the expected opposite-strand fraction
under gDNA equals the expected same-strand fraction — the term becomes
uninformative. This is the correct behaviour: the pathway
self-silences when it has no discriminating power.

---

## 5. Unified joint-likelihood calibration

### 5.1. The objective

Combine all three information sources into a single log-likelihood:

$$\ell(\lambda_G) = \underbrace{\ell_\text{density}(\lambda_G)}_\text{§2.1: Poisson mixture over regions} \;+\; \underbrace{\ell_\text{strand}(\lambda_G)}_\text{§4: opposite-strand exonic} \;+\; \underbrace{\ell_\text{unmapped}(\lambda_G)}_\text{§3.4: uT:3 Poisson}$$

Maximise via 1-D Newton–Raphson on $\log\lambda_G$ (log-concave,
typically converges in 3–5 iterations).

### 5.2. Properties

1. **No arbitrary cutoffs.** No P₁₀, no blend weights, no hard
   mappability mask.
2. **Self-weighting by Fisher information.** Each pathway contributes
   proportional to its effective exposure (span) and signal strength.
   In stranded libraries, the strand pathway dominates. In unstranded
   libraries, the unmapped pathway dominates. The density pathway
   fills in the middle.
3. **Graceful degradation.** If `uT` tags are absent (e.g.
   non-STAR aligner), the unmapped term drops out. If no mappable BED
   is provided, $f_i = 1$ for all regions. If $\sigma \approx 0.5$,
   the strand term contributes nearly nothing. Each pathway failing is
   equivalent to removing a term from the likelihood — the remaining
   terms still produce a valid estimate.
4. **Rich diagnostics.** Per-pathway λ_G estimates, Fisher information
   shares, inter-pathway consistency (chi-squared test on the three
   point estimates) — all reportable in `summary.json`.
5. **The zero-count moment estimator (§2.2)** serves as an
   always-valid lower bound and consistency check.

### 5.3. Expected behaviour on the pipeline-validation BAM

On the unstranded + 50% gDNA STAR run from April 18:

| Pathway | Available signal | Expected contribution |
|---|---|---|
| Strand | σ = 0.500 → uninformative | ~0% of Fisher info |
| Density (mixture MLE) | 683k regions, ~54% zeros | Modest; constrains λ_G from above via zeros |
| Unmapped (uT:3) | ~300k–1.5M reads over ~300 Mb | **Dominant**; pins λ_G precisely |

The unmapped pathway is predicted to rescue calibration in exactly the
regime where the current system fails.

---

## 6. Implementation plan

### Phase 1: Data collection (unmapped pathway prototype)

**Goal:** Validate that uT:3 counts produce a credible λ_G on the
existing pipeline-validation BAM.

1. Write a diagnostic script (`scripts/debug/unmapped_pathway_probe.py`)
   that:
   - Reads the name-sorted STAR BAM via pysam
   - Counts fragments by `uT` tag value (0–4)
   - Computes $S_\text{rep}$ from
     `alignable_unmappable_rl125_grch38_star.bed`
   - Estimates $\hat\lambda_G = (1-\pi_\text{RNA})\,N_\text{uT3} / S_\text{rep}$
   - Compares to truth and to the existing pathway estimates
2. Run on the April 18 BAM:
   `/scratch/.../mappability_stress/runs/gdna_high_ss_0.50_nrna_none/star/name_sorted.bam`
3. If λ_G is within 20% of truth → proceed to Phase 2.

**Estimated effort:** 1–2 hours.

### Phase 2: Poisson mixture density estimator

**Goal:** Replace P₁₀ with a mixture-MLE density pathway.

1. Add `calibration_mixture.py` (or extend `calibration.py`) with:
   - `poisson_mixture_mle(counts, spans, mappable_fracs, max_iter=50)`
     → returns `(lambda_g, pi_rna_per_region, fisher_info)`
   - `zero_count_bound(counts, spans, mappable_fracs)` → returns
     `lambda_g_lower_bound`
2. Compute $f_i$ (mappable fraction) at index time:
   - In `index.py`, for each region, compute fraction of span covered
     by mappable intervals → store as a column in `regions.feather`
   - This replaces the binary `compute_containment_mask()` with a
     continuous weight
3. Wire into `calibration.py`:
   - Replace `_density_pathway()` internals with the mixture MLE
   - Keep the function signature compatible
4. Unit tests:
   - Synthetic Poisson data with known λ_G → verify recovery
   - Edge cases: all zeros, single region, $f_i = 0$ regions
5. Run existing test suite (`pytest tests/ -v`) — verify no regression.

**Estimated effort:** 1–2 days.

### Phase 3: C++ scanner changes (uT tag counting)

**Goal:** Make unmapped-read counts available to the Python calibration
layer.

1. In `bam_scanner.cpp`, add counters for each `uT` tag value (0–4)
   to the scan result struct. For unmapped reads (BAM flag 0x4), read
   the `uT` aux tag and increment the appropriate counter.
2. Expose via nanobind in `_bam_impl` → `ScanResult.unmapped_counts`
   (dict or array of 5 ints).
3. In `buffer.py` or `scan.py`, propagate the counts to the Python
   layer.
4. Store $S_\text{rep}$ (total unmappable span) in the index metadata
   at build time (from the same BED file used for mappable regions).

**Estimated effort:** half a day (small C++ change, mostly plumbing).

### Phase 4: Joint-likelihood unification

**Goal:** Replace the two-pathway blend with the unified objective.

1. New function `calibrate_gdna_joint()` in `calibration.py`:
   - Inputs: region counts/spans/mappable-fracs (density), opp-strand
     counts + exonic span + σ (strand), uT:3 count + $S_\text{rep}$
     + $\pi_\text{RNA}$ (unmapped)
   - 1-D Newton on $\log\lambda_G$
   - Returns: `CalibrationResult` with per-pathway diagnostics
2. Deprecate (but keep) `calibrate_gdna()` for backward compatibility.
3. Wire into `pipeline.py` as the default calibration path.
4. Update `summary.json` schema:
   - `calibration.lambda_gdna` (joint MLE)
   - `calibration.density_pathway_lambda` (mixture MLE alone)
   - `calibration.strand_pathway_lambda` (strand alone)
   - `calibration.unmapped_pathway_lambda` (unmapped alone)
   - `calibration.unmapped_uT_counts` (array of 5)
   - `calibration.zero_count_lambda_lower_bound`
   - `calibration.fisher_info_shares` (3-element: density, strand,
     unmapped)
5. Full test suite + benchmark sweep.

**Estimated effort:** 1–2 days.

### Phase 5: Validation

**Goal:** Confirm improvement across the design space.

1. Re-run the 24-scenario analytical stress test with the new
   calibration (requires extending the synthetic framework to generate
   uT:3-like counts).
2. Re-run the STAR pipeline validation at SS = {0.5, 0.9, 1.0} ×
   gDNA = {low, high}.
3. Re-run the full-scale hulkrna benchmarks
   (`scripts/benchmarking/`).
4. Compare λ_G accuracy, mRNA/gDNA/nRNA totals, and per-gene
   correlation against truth.
5. Publish results to `docs/mappability/unified_validation.md`.

**Estimated effort:** 1–2 days (mostly compute time).

---

## 7. Risk assessment

| Risk | Severity | Mitigation |
|---|---|---|
| π_RNA prior is wrong for some tissues (e.g. testes with high LINE-1 RNA) | Medium | Make π_RNA a CLI parameter; report per-pathway λ_G so discrepancies are visible; use Beta prior with proper uncertainty propagation |
| Non-STAR aligners don't emit `uT` tags | Low | Unmapped pathway gracefully drops out; density + strand pathways still work. Minimap2 doesn't produce `uT` but does report unmapped reads with flag 0x4 — could fall back to total-unmapped with a wider π_RNA |
| Poisson mixture EM doesn't converge on mega-loci with extreme RNA dominance | Low | Profile likelihood is log-concave in λ_G; convergence is guaranteed. Initialise from the zero-count bound |
| $S_\text{rep}$ depends on read length | Low | Variation is <5% between 100–150 bp reads. Store a single value; optionally allow user override |
| Backward compatibility of `summary.json` | Low | New fields are additive; old fields retained |

---

## 8. Open questions for discussion

1. **Should the mixture MLE use region-level or equivalence-class-level
   counts?** Region-level is simpler and matches the current
   architecture. EC-level would be more principled but requires deeper
   plumbing changes.

2. **Should $f_i$ be computed from the mappable BED or from observed
   multimapper density?** The BED is a static property of the genome +
   aligner; the observed multimapper density is data-adaptive but
   potentially circular. Starting with the BED is safer.

3. **Should we attempt to decompose uT:2 (too many mismatches)?**
   Probably not — the contamination signal is too mixed. Better to
   report it as a QC diagnostic and flag libraries with uT:2 > 3% as
   potentially contaminated.

4. **Is there value in a second-pass that re-examines uT:3 reads with
   relaxed STAR settings?** STAR's `--outSAMunmapped Within` can
   report placeholder coordinates for unmapped reads. This could
   let us see if uT:3 reads cluster in repeat regions (gDNA) vs.
   exonic repeats (RNA). This is an enhancement, not a prerequisite.

5. **Can we calibrate π_RNA from stranded libraries?** In a stranded
   library, repeat-element RNA should show strand bias while gDNA
   should not. A training set of stranded libraries could provide a
   tissue-specific π_RNA lookup table. Worth investigating but not
   blocking.

---

## 9. Summary

| Problem | Current approach | Failure mode | Proposed replacement |
|---|---|---|---|
| Estimating λ_G from region densities | P₁₀ percentile | Collapses to zero under Poisson sparsity | Poisson mixture MLE with continuous mappable weights |
| Combining strand and density signals | Heuristic blend weight w_strand | Magic constant; silent failure when one pathway is uninformative | Joint log-likelihood; self-weighting by Fisher information |
| Unmapped reads | Ignored | 3–15% of library discarded; dominant gDNA signal lost | uT:3 Poisson pathway with unmappable-genome denominator |
| Mappability mask | Hard in/out containment filter | Loses 49% of regions; can hurt when aligner retains most reads | Continuous mappable fraction $f_i \in [0,1]$ per region |

The core insight: **every source of gDNA information becomes a term in
the same likelihood, weighted by its natural Fisher information, with
no magic constants.** The P₁₀ goes away, the blend weight goes away,
the hard mask goes away, and the unmapped reads stop being invisible.
