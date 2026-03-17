# gDNA Calibration: Expression-Based Purity Scoring Plan

## Motivation

The calibration module estimates global gDNA nuisance parameters (strand
overdispersion κ, fragment-length distribution, background density) from
region-level evidence.  The current seed selection uses hard thresholds
(zero spliced, strand symmetry z-test, low-density percentile), which
fail under sparse simulation conditions and conflate density with purity.

### Critical Constraint: Cannot Rely on Intergenic Regions Alone

Many RNA-seq libraries use hybrid capture (exome capture) which **depletes
intergenic reads**.  Capture probes pull down DNA and RNA indiscriminately,
producing high exonic density even for gDNA.  Intergenic regions CAN be
included but MUST NOT be the only seed source.

## Core Principle

Every region partition region — intergenic, intronic, or exonic — gets
scored for **RNA expression evidence**.  Regions with the least evidence
are the purest gDNA candidates.  This works for:

- **Standard polyA / rRNA-depletion libraries**: Intergenic and
  unexpressed-gene regions both qualify as seeds.
- **Hybrid capture libraries**: Intergenic depleted, but unexpressed
  exonic/intronic regions within capture space still carry gDNA signal.
- **Any mixture**: The purity score naturally adapts.

## Algorithm Phases

### Phase 1 — Compute Region Stats (expand existing)

Add annotation flags to the stats dict:

```python
stats["tx_pos"]   = region_df["tx_pos"].values.astype(bool)
stats["tx_neg"]   = region_df["tx_neg"].values.astype(bool)
stats["exon_pos"] = region_df["exon_pos"].values.astype(bool)
stats["exon_neg"] = region_df["exon_neg"].values.astype(bool)
stats["ref"]      = region_df["ref"].values
```

### Phase 2 — Compute RNA Purity Score (new)

Each region receives a composite penalty score.  **Higher = more evidence
of RNA expression** (i.e. less pure as gDNA).  Uses three additive
signals, all computable without knowing gDNA density:

| Signal | Measures | Method |
|--------|----------|--------|
| **Splice penalty** | Spliced reads = definitive RNA | `splice_rate` scaled by confidence |
| **Strand penalty** | Sense excess over symmetric | Binomial z-score for p ≠ 0.5, weighted by SS |
| **Density penalty** | Coverage excess over background | z-score vs median density of candidate regions |

The density background for the z-score is bootstrapped from all
zero-spliced, approximately symmetric regions (loose pre-filter).
This avoids circular dependency — no gDNA density estimate needed.

Regions are **sorted by purity score**; the bottom percentile (purest)
become seeds.  All region types (intergenic, intronic, exonic) are
eligible.

### Phase 3 — Estimate Initial Parameters from Seed

- **gDNA density**: `Σ n_unspliced[seed] / Σ region_length[seed]`.
  **FROZEN** — not updated in the iteration loop.
- **κ_sym**: MoM from seed strand ratios (existing `estimate_kappa_sym`).
- **FL model**: FL observations from seed regions only.

### Phase 4 — Iterate (κ and weights only)

```
for iteration in 1..max_iter:
    weights = score_regions(stats, kappa_sym, SS, frozen_density)
    new_kappa = estimate_kappa_sym(stats, weights)
    if converged: break
    kappa_sym = new_kappa
```

Density is frozen.  Only κ and per-region weights evolve.  This breaks
the circular dependency that caused κ to hit bounds at SS = 0.5.

### Phase 5 — Final FL Model

Build gDNA FL from FL observations in regions where `weight > 0.8`.
With frozen density anchoring the weights, high-weight regions are
genuinely gDNA-enriched.

### Phase 6 — Early Bailout

If fewer than `MIN_SEED_REGIONS` (e.g. 5) pass the minimum purity
threshold, return a "no gDNA detected" calibration with density = 0.

## Hierarchical Extension (future)

Wire `stats["ref"]` for per-reference grouping:

1. Per-ref density from seed regions on that ref.
2. Shrink toward global: `ref_density = w × raw + (1−w) × global`.
3. Feed per-ref density into `score_regions()` as a per-region array.

Mirrors `_compute_ref_gdna_densities` in `locus.py`.

## Code Changes

| Component | Status | Change |
|-----------|--------|--------|
| `compute_region_stats()` | Expand | Add annotation & ref arrays |
| `select_seed()` | **Replace** | Purity-score-based ranking |
| `estimate_kappa_sym()` | Keep | Works with weighted regions |
| `estimate_gdna_density()` | Keep | Used once from seed, frozen |
| `score_regions()` | Keep | Receives frozen density |
| `build_gdna_fl_model()` | Keep | Fed high-weight regions |
| `calibrate_gdna()` | Modify | Freeze density, bailout, new seed |

## Simulation Improvements

1. **Scale up**: 5M genome, 200+ transcripts, 100k fragments.
2. **Hybrid capture scenario**: Deplete intergenic reads to validate
   purity score finds seeds in unexpressed exonic regions.
3. **Vary RNG seed** across runs for statistical power.
