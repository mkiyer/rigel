# Plan: Anchored Empirical Bayes for gDNA Resolution

**Status**: Proposal  
**Date**: 2026-04-07  
**Prerequisite**: Dirichlet-Multinomial FL prior (implemented), fraction-based priors (implemented)

## Context

At SS=0.50, gDNA and nRNA are degenerate in likelihood space — identical strand signal, near-identical FL, and (for single-transcript loci) identical effective lengths. The current prior budget `C=1.0` is too weak to break the tie. This plan addresses both sides of Bayes' equation: **Prior** (Idea 3) and **Likelihood** (Idea 4).

---

## Part A: Informative Prior (Idea 3)

### Current State

```python
# locus.py compute_locus_priors()
gamma = e_sum / n_sum     # local gDNA mixing fraction from calibration
alpha_gdna = gamma * C    # C = 1.0
alpha_rna  = (1 - gamma) * C
```

With `C=1.0` and 2000 fragments, the prior has 0.05% influence. It's noise.

### Proposed Change: Scale C by Locus Fragment Count

The simplest principled approach: make the prior *competitive* with the data by scaling C proportionally to the number of fragments entering the locus EM:

$$C_{\text{locus}} = \kappa \cdot N_{\text{locus}}$$

where $\kappa \in [0, 1]$ controls the prior's relative weight. This gives:

$$\alpha_{\text{gDNA}} = \gamma \cdot \kappa \cdot N_{\text{locus}}$$
$$\alpha_{\text{RNA}} = (1 - \gamma) \cdot \kappa \cdot N_{\text{locus}}$$

#### Choosing κ: Strand Specificity as a Natural Weight

The degree of strand specificity directly controls how much we trust the calibration's gDNA estimate:

- **SS=1.0**: Strand decomposition is exact → calibration is highly reliable → strong prior
- **SS=0.50**: Only the density baseline is available → calibration is noisy → moderate prior (but still needed because the likelihood is degenerate!)

The existing blending weight `w = (2·SS − 1)²` already captures this trust gradient. However, at SS=0.50, `w=0` — and this is precisely where we *need* the prior most. The relationship should actually be **inverted**: the prior should be strongest when the likelihood is weakest.

**Proposed formula:**

$$\kappa = \kappa_{\min} + (1 - w) \cdot (\kappa_{\max} - \kappa_{\min})$$

where:
- `w = (2·SS − 1)²` (strand reliability weight)
- `κ_min`: Minimum prior strength for perfectly stranded data (e.g., 0.01 — the calibration is good but likelihood can override)
- `κ_max`: Maximum prior strength for unstranded data (e.g., 0.5 — prior locks in roughly half the evidence weight)

At SS=0.50: `w=0`, `κ = κ_max = 0.5`. With N=2000 and γ=0.01: `α_gDNA = 0.01 × 0.5 × 2000 = 10`, `α_RNA = 0.99 × 0.5 × 2000 = 990`. The prior "expects" 10 gDNA fragments out of 2000, with the strength of ~1000 pseudo-observations. The EM can still deviate if the evidence is strong, but for the degenerate case (identical likelihoods), the prior holds.

At SS=1.0: `w=1`, `κ = κ_min = 0.01`. Prior = 20 pseudo-observations on 2000 fragments. Light touch — the strand signal in the likelihood already separates gDNA from RNA.

#### Alternative: Direct Expected Count as Prior

An even simpler approach avoids the fraction-based formula entirely. From calibration we have `λ_G` (gDNA density in frags/bp) and `locus_span`. The expected gDNA count is:

$$E[N_{\text{gDNA}}] = \lambda_G \cdot L_{\text{locus}}$$

We could use this directly as the Dirichlet pseudocount:

$$\alpha_{\text{gDNA}} = \kappa \cdot \lambda_G \cdot L_{\text{locus}}$$
$$\alpha_{\text{RNA}} = \kappa \cdot (N_{\text{locus}} - \lambda_G \cdot L_{\text{locus}})$$

This is mathematically equivalent to the γ-fraction approach (since γ ≈ λ_G · L / N), but is more interpretable: the prior says "expect this many gDNA fragments," with strength scaled by κ.

#### Advantages

- **Simple**: One new parameter (κ), or two (κ_min, κ_max) for the SS-dependent version.
- **Uses existing architecture**: Only `compute_locus_priors()` changes. No EM restructuring.
- **Principled**: κ scales with 1−w, so the prior is strongest when the likelihood is weakest.
- **Data-driven**: γ comes from calibration (thousands of regions, millions of fragments).

#### Risks

- **γ miscalibration propagates**: If calibration over- or under-estimates γ, the error is amplified by N_locus. Mitigation: κ < 1.0 (never let the prior completely dominate).
- **No local gDNA variation**: The prior forces all loci toward the global γ. True local spikes (e.g., amplified regions) will be attributed to RNA. For the unstranded degenerate case, this is acceptable — there's no way to distinguish local gDNA spikes from high nRNA without strand signal.

---

## Part B: Likelihood Differentiation (Idea 4)

### Sub-idea 4a: Remove Fragment Length Subtraction for gDNA (Simplest)

**Current gDNA bias correction:**
```
eff_len = max(locus_span - genomic_footprint + 1, 1)
log_lik -= log(eff_len)
```

This assumes gDNA fragments are drawn uniformly from a region of size `locus_span`, just like nRNA fragments are drawn from their genomic span. For a single-transcript locus with a 10kb genomic span and 200bp fragments:

- nRNA eff_len = 10000 − 200 + 1 = 9801
- gDNA eff_len = 10000 − 200 + 1 = 9801 → **identical**

**Proposed**: For gDNA, use the raw span without fragment length subtraction:

```
gDNA eff_len = locus_span    (no subtraction)
```

**Biological rationale**: gDNA fragments represent sheared genomic DNA. The Library preparation process fragments DNA randomly, then fragments are captured by size selection. The "effective length" concept (positions where a fragment of length F can start) applies to RNA because transcripts have defined boundaries. But gDNA is continuous — a fragment at position X could equally well have been a fragment at position X+1 that was sheared slightly differently. The gDNA component represents the *density* of a uniform Poisson process, not a transcript with boundaries.

**Impact at the numbers**: With locus_span=10000 and frag=200:

- nRNA: `−log(9801) = −9.19`
- gDNA: `−log(10000) = −9.21`
- Difference: **−0.02 log-units per fragment**

This is small. Over 2000 fragments: `2000 × 0.02 = 40` log-units total — enough to shift the EM but only marginal per-fragment discrimination.

**Verdict**: Correct in principle but insufficient alone. The fragment-length subtraction is a second-order effect (~2% of log).

### Sub-idea 4b: Extend gDNA Span Beyond Locus Boundaries

**Current**: `bias_profiles[gdna_idx] = locus_span` (union of transcript genomic spans within the locus).

**Proposed**: Extend the gDNA span to include flanking intergenic sequence:

```
gDNA_span = locus_span + 2 × flank_size
```

**Biological rationale**: gDNA truly is uniform across the genome. A locus at chr1:1000000-1010000 will have gDNA fragments at chr1:999800 and chr1:1010200 — they just don't overlap any transcript and are currently discarded as intergenic. The gDNA component's *true* effective length extends indefinitely, but practically, fragments that don't overlap the locus can't compete in the EM. Still, the gDNA effective length should reflect that the gDNA density extends beyond the locus boundaries.

**Impact**: With 10kb locus + 5kb flanking on each side:

- nRNA span: 10kb → `−log(9801) = −9.19`
- gDNA span: 20kb → `−log(19801) = −9.89` (or `−log(20000) = −9.90` without F subtraction)
- Difference: **−0.70 log-units per fragment**

That's 35× more discrimination than Sub-idea 4a. Over 2000 fragments, this is a 1400 log-unit total shift — decisive.

**Choice of flank size**: The natural scale is the fragment length distribution itself. At the tails of the FL distribution (e.g., 500bp), a fragment centered 250bp outside the locus boundary would still partially overlap and be "captured." Using `flank = max_fragment_length / 2` is defensible. Alternatively, a fixed value like 1000bp or 5000bp could work.

**But there's a cleaner formulation**:

### Sub-idea 4c: Use Global gDNA Effective Length (Chromosome or Genome Scale)

**The "Nuclear Option"**: gDNA fragments are drawn from the *entire genome*. The gDNA component's effective length is `G` (genome size, ~3 Gbp for human). This gives:

```
gDNA eff_len = genome_size    (~3e9)
```

- nRNA: `−log(9801) = −9.19`
- gDNA: `−log(3e9) = −21.82`
- Difference: **−12.63 log-units per fragment** → massive penalty

This is too strong — it would effectively zero out gDNA everywhere. The gDNA component could never compete.

**Why this fails**: The gDNA likelihood already embeds the genome-scale normalization implicitly. The gDNA density `λ_G` is defined as frags/bp, already normalized by genome size. Adding `−log(G)` double-counts the normalization.

**Verdict**: Reject 4c. The genome-wide normalization belongs in the prior (`α_gDNA = λ_G × span`), not in the per-fragment likelihood.

### Sub-idea 4d: Incorporate Fragment Overhang Information

**Observation**: The C++ scorer computes `tx_end = tx_start + frag_len` for RNA components, which can exceed the transcript length (overhang). For gDNA, `tx_start=0` and `tx_end=footprint`. Fragments that overhang the transcript boundary are more likely to be gDNA or nRNA (arising from the unspliced pre-mRNA or genomic DNA that extends beyond the annotated transcript).

**Current behavior**: Overhang is penalized by the `oh_log_penalty` (default `log(0.01) ≈ −4.6` per overhanging base). This already discriminates against mRNA for overhanging fragments. However, nRNA components (synthetic single-exon transcripts spanning the full genomic region) typically don't overhang — the fragment falls within the nRNA span.

**What about fragments that overhang the nRNA span itself?** These would be fragments extending *beyond* the TSS or TES of the merged transcript group. Such fragments:
- Can't be nRNA (nRNA is bounded by TSS/TES)
- Can't be mRNA (mRNA is bounded by TSS/TES exons)
- **Must be gDNA** (unconstrained)

Currently: these fragments get scored against both nRNA and gDNA with similar likelihoods (the overhang penalty applies to nRNA, the gDNA component absorbs them). The issue is that in the degenerate case, fragments *within* the span can't tell gDNA from nRNA.

**Verdict**: Overhang is already handled correctly (oh_penalty discriminates). The degenerate problem is fragments *within* the common span, not at boundaries.

### Sub-idea 4e: Absorb Flanking Intergenic Fragments into Locus gDNA

**Current**: Intergenic fragments (no transcript overlap) are counted in pipeline statistics but never enter the EM. They are effectively "wasted" — assigned to gDNA in the summary but without influencing the EM's gDNA estimate.

**Proposed**: For each locus, extend the fragment capture zone by some flanking distance. Intergenic fragments within this zone are added to the locus EM as gDNA-only candidates (they can only be gDNA since they don't overlap any transcript).

**Impact**:
1. These fragments provide *pure gDNA signal* — no ambiguity with mRNA or nRNA
2. They directly augment the gDNA component's EM count, anchoring it to a physical reality
3. The number of flanking intergenic fragments is proportional to `λ_G × flank_size × 2` — exactly the gDNA density information we want

**Advantage**: No parameter tuning for effective length. The data itself tells the EM how much gDNA there is, calibrated by the physical density of intergenic fragments near the locus.

**Disadvantage**: Requires buffering and routing intergenic fragments to nearby loci, which adds complexity to the scan/scoring pipeline. The C++ resolve step currently skips intergenic fragments entirely — they'd need to be assigned to the nearest locus.

### Sub-idea 4f: Merge Spatially Adjacent Loci with Shared gDNA

**Concept**: Instead of each locus having its own independent gDNA component, merge physically adjacent loci (within some distance threshold) into a single EM problem with one shared gDNA component. The shared gDNA component's effective length spans the entire merged region (including the intergenic gaps between loci).

**Advantage**: The shared gDNA density is constrained across loci. If locus A has high unspliced density and neighboring locus B has low density, the shared gDNA component can't simultaneously explain both — the excess in A must be nRNA.

**Disadvantage**: This is essentially building "mega-loci" artificially, which we already know cause performance issues. The EM complexity grows with the number of components. Also, the motivation for merging is specifically the gDNA degeneracy — at SS>0.7, independent loci work fine. This would add complexity for a problem that only manifests in unstranded data.

---

## Evaluation Matrix

| Sub-idea | Discrimination per fragment | Complexity | Biological soundness | Risk |
|----------|---------------------------|------------|---------------------|------|
| **4a**: No F subtraction for gDNA | −0.02 log-units | Trivial | Defensible | Low |
| **4b**: Extend span with flanking | −0.70 log-units (5kb flank) | Low | Strong | Flank size choice |
| **4e**: Absorb intergenic frags | Data-driven (exact) | Medium-high | Strongest | Pipeline complexity |
| **4f**: Merge adjacent loci | N/A (structural) | High | Moderate | Performance, over-engineering |
| **4c**: Genome-scale eff_len | −12.6 log-units | Trivial | Wrong (double-counts) | **Reject** |
| **4d**: Overhang info | Already handled | N/A | N/A | N/A |

---

## Recommended Plan

### Phase 1: Informative Prior (Idea 3) — Immediate

**Change `compute_locus_priors()` to scale C by locus fragment count:**

```python
def compute_locus_priors(
    loci, index, calibration,
    n_locus_fragments,           # NEW: array of per-locus fragment counts
    kappa_min=0.01,              # NEW: prior strength at SS=1.0
    kappa_max=0.5,               # NEW: prior strength at SS=0.50
):
    w = (2 * SS - 1) ** 2
    kappa = kappa_min + (1 - w) * (kappa_max - kappa_min)

    for li, locus in enumerate(loci):
        gamma = ...  # unchanged
        C_locus = kappa * n_locus_fragments[li]
        alpha_gdna[li] = gamma * C_locus
        alpha_rna[li]  = (1 - gamma) * C_locus
```

**Scope**: Python-only change to `locus.py` and `pipeline.py` (pass fragment counts). No C++ changes. No EM restructuring.

**Testing**: Existing scenario tests at SS=0.65 and SS=1.0. Add new unstranded (SS=0.50) scenario.

### Phase 2: gDNA Effective Length (Idea 4a + 4b) — Immediate

**Two changes, both in `locus.py`:**

1. **4a**: Remove fragment length subtraction for gDNA by setting `tx_start = tx_end = 0` instead of `tx_start=0, tx_end=footprint`. This makes `frag_len=0` in the bias correction, so `eff_len = locus_span` (the raw span).

2. **4b**: Extend `gdna_span` by a configurable flanking distance. In `build_loci()`, when computing `gdna_span`, add flanking:

```python
# After computing merged_intervals and span:
gdna_span = max(span + 2 * flank, 1)
```

The flank distance should be configurable (default: the 99th percentile of the global FL distribution, which is a data-driven proxy for "how far outside the locus can gDNA fragments reach").

**Note on 4a**: There is actually a simpler way to achieve this. Instead of changing the CSR data (`tx_start`/`tx_end`), change `bias_profiles[gdna_idx]` to use `gdna_span + flanking` *without* changing the footprint subtraction. The effect is:

```
eff_len = (gdna_span + flank) - footprint + 1
```

versus nRNA:
```
eff_len = nrna_span - footprint + 1
```

This naturally creates a differential of `flank` in the denominator. For 5kb flank on a 10kb locus: gDNA eff_len = 15801 vs nRNA eff_len = 9801. Per fragment: `−log(15801) + log(9801) = −0.48` log-units advantage to nRNA.

**Scope**: Python-only changes to `locus.py`. Possibly a new config parameter for flank distance.

### Phase 3: Intergenic Fragment Utilization (Idea 4e) — Future

**Deferred** because it requires changes to the C++ resolve/scoring pipeline and fragment routing. This would be a meaningful improvement but is architecturally more invasive. The combination of Phase 1 + Phase 2 should be sufficient for the immediate problem.

If Phase 1+2 are insufficient, this becomes the next priority. The implementation would:
1. In the C++ scanner, buffer intergenic fragments with their genomic position
2. After locus construction, assign each intergenic fragment to the nearest locus within a configurable distance
3. Add these as gDNA-only candidates in the CSR data

### Phase 4: Locoregional Density Prior — Future

Instead of using the global γ, compute a locoregional density from a wider neighborhood (e.g., 1Mb window) using the region partition. This would capture copy-number-driven local variation in gDNA density. The region partition already provides per-region E[gDNA] — summing over a sliding window of regions is straightforward.

**Deferred** because the global γ is a good starting point, and local variation in gDNA density is a second-order effect compared to the gDNA/nRNA degeneracy.

---

## Parameter Inventory

| Parameter | Location | Default | Purpose |
|-----------|----------|---------|---------|
| `kappa_min` | `CalibrationConfig` | 0.01 | Prior strength at SS=1.0 (light touch) |
| `kappa_max` | `CalibrationConfig` | 0.50 | Prior strength at SS=0.50 (strong anchor) |
| `gdna_flank` | `CalibrationConfig` | `fl_p99` or 1000 | Flanking extension for gDNA effective length (bp) |

These replace the current `total_pseudocount = 1.0`.

---

## Summary

**Phase 1** (Prior): Makes the calibrated gDNA fraction a real constraint by scaling the prior budget with locus size and inversely with strand specificity. Simple, immediate impact.

**Phase 2** (Likelihood): Creates per-fragment discrimination by extending gDNA's effective length beyond the locus boundary. Biologically motivated — gDNA truly extends beyond transcript boundaries. Complements the prior by giving the EM actual likelihood signal.

Together, these attack the degeneracy from both sides of Bayes' theorem. The prior says "expect this much gDNA." The likelihood says "each fragment is slightly more probable under nRNA than gDNA." In the degenerate case (SS=0.50), the prior dominates — which is correct, because the global density calibration is the best information available. In the stranded case (SS=1.0), the likelihood dominates — the prior is just a light touch.
