# Unified EM Initialization: Density + Strand Hybrid Model

## Status: Phase 1a/1b/1d IMPLEMENTED — Phase 1c pending

## 1. Problem Statement

The rigel EM solver must initialize three pools per locus: **mRNA**
(mature RNA), **nRNA** (nascent/pre-mRNA), and **gDNA** (genomic DNA).
The current initialization has three critical flaws:

1. **Hard strand threshold**: The pipeline requires strand specificity
   ≥ 0.6 and refuses to run unstranded data, excluding a large fraction
   of public RNA-seq datasets.

2. **gDNA-blind nrna_frac prior**: The nascent fraction prior (`nrna_frac_t`) uses raw
   unspliced/spliced exonic ratios without subtracting gDNA background,
   inflating nrna_frac for silent genes and single-exon transcripts.

3. **No density information**: Coverage density differences between exonic,
   intronic, and intergenic space are not used for initialization, despite
   being a powerful gDNA-independent signal available in ALL libraries.

This document specifies a **unified hybrid initialization model** that
combines strand-based and density-based evidence through inverse-variance
weighting, enabling rigel to handle the full spectrum from perfectly
stranded to completely unstranded libraries.

## 2. Theoretical Foundation

### 2.1 The Three Spatial Zones

Genomic space is partitioned into three zones with distinct fragment
compositions:

| Zone | Contains | Density |
|------|----------|---------|
| **Intergenic** | gDNA only | D_ig = λ_gDNA |
| **Intronic** | gDNA + nRNA | D_in = λ_gDNA + λ_nRNA |
| **Exonic** | gDNA + nRNA + mRNA | D_ex = λ_gDNA + λ_nRNA + λ_mRNA |

where λ denotes read density (fragments per base pair).

### 2.2 Density Subtraction (unstranded signal)

By subtracting adjacent zone densities, biological rates are isolated:

```
λ_gDNA = D_intergenic
λ_nRNA = max(0, D_intronic − λ_gDNA)
λ_mRNA = max(0, D_exonic − D_intronic)
```

**Key property**: gDNA cancels in the mRNA estimate (D_exonic − D_intronic
removes both gDNA and nRNA background).

**Limitation**: Assumes uniform gDNA density across zones. Capture-enriched
libraries violate this (exonic gDNA density ≫ intronic gDNA density).

### 2.3 Strand Subtraction (stranded signal)

For a library with strand specificity s ∈ [0.5, 1.0]:

- RNA produces s sense + (1−s) antisense fragments
- gDNA produces 0.5 sense + 0.5 antisense (double-stranded, no strand bias)

Observed counts in any zone:
```
Sense     = R × s + G × 0.5
Antisense = R × (1−s) + G × 0.5
```

Subtracting:
```
Sense − Antisense = R × (2s − 1)
R = (Sense − Antisense) / (2s − 1)
```

**Key property**: gDNA cancels completely regardless of absolute level.
Works even with capture-enriched exonic gDNA.

**Limitation**: Denominator (2s−1) → 0 as s → 0.5 (unstranded), causing
the estimate to explode.

### 2.4 Unified Hybrid Estimator

The two estimators have complementary failure modes:
- **Density**: works for unstranded, fails for capture-enriched
- **Strand**: works for capture-enriched, fails for unstranded

We combine them via **inverse-variance weighting**. The variance of the
strand estimator is proportional to 1/(2s−1)², so its reliability weight
is:

```
W_strand = (2s − 1)²     ∈ [0, 1]
W_density = 1 − W_strand
```

| Library type | s | W_strand | W_density | Behavior |
|---|---|---|---|---|
| Perfect stranded | 1.0 | 1.0 | 0.0 | Pure strand subtraction |
| Good stranded | 0.95 | 0.81 | 0.19 | Mostly strand |
| Weak stranded | 0.7 | 0.16 | 0.84 | Mostly density |
| Unstranded | 0.5 | 0.0 | 1.0 | Pure density subtraction |

The unified rate estimates for any zone are:

```
λ_final = W_strand × λ_strand + W_density × λ_density
```

## 3. Initialization Quantities

The EM requires initialization of three sets of values:

### 3.1 Per-locus gDNA (γ_ℓ)

**Already well-implemented** via Empirical Bayes shrinkage of intergenic
density. The current `compute_eb_gdna_priors()` estimates a global gDNA
rate from intergenic fragments, then shrinks per-locus estimates toward
this global rate. This remains unchanged.

### 3.2 Per-transcript θ_t (total RNA abundance)

Currently initialized from unambig exonic counts. Under the linked model,
θ_t = mRNA_t + nRNA_t. The density/strand subtraction improves this by
providing a gDNA-corrected estimate of total RNA.

### 3.3 Per-transcript nrna_frac_t (nascent RNA fraction)

This is where the hybrid estimator is essential. The nrna_frac prior cascade
(transcript → TSS group → locus-strand) uses the unified hybrid
estimator at each level.

## 4. Detailed Algorithm

### 4.1 Inputs (available at initialization time)

From the scan phase, we have per-transcript unambig fragment counts broken
down by splice status and strand:

```
For each transcript t:
  spliced_sense[t]       — spliced exonic, sense strand
  spliced_anti[t]        — spliced exonic, antisense strand
  unspliced_sense[t]     — unspliced exonic, sense strand
  unspliced_anti[t]      — unspliced exonic, antisense strand
  intronic_sense[t]      — intronic, sense strand
  intronic_anti[t]       — intronic, antisense strand
```

From the index:
```
  L_exonic[t]            — total exonic length (bp)
  L_intronic[t]          — total intronic span (bp), 0 for single-exon
```

From the gDNA EB prior:
```
  gdna_density           — global intergenic density (frags / bp)
  locus_gdna_density[ℓ]  — per-locus EB-shrunk gDNA density
```

From the strand model:
```
  s                      — strand specificity (p_sense for the exonic
                           spliced sub-model, range [0.5, 1.0])
```

### 4.2 Hybrid Rate Estimation (per transcript)

For each transcript t, compute mRNA and nRNA rate estimates:

```python
def compute_hybrid_rates(t, s, gdna_density):
    """Returns (lambda_mrna, lambda_nrna, evidence_weight)."""

    # ── Strand weight ──
    denom = 2.0 * s - 1.0
    W = denom ** 2 if denom > 0.01 else 0.0

    # ── Total counts ──
    exon_total = (spliced_sense[t] + spliced_anti[t]
                  + unspliced_sense[t] + unspliced_anti[t])
    intron_total = intronic_sense[t] + intronic_anti[t]

    # ── Density estimator ──
    if L_exonic[t] > 0:
        D_exon = exon_total / L_exonic[t]
    else:
        D_exon = 0.0

    if L_intronic[t] > 0:
        D_intron = intron_total / L_intronic[t]
    else:
        # Single-exon transcript: no intronic signal.
        # nRNA rate is unobservable from density alone.
        D_intron = gdna_density  # assume intron = background

    den_nrna = max(0.0, D_intron - gdna_density)
    den_mrna = max(0.0, D_exon - D_intron)

    # ── Strand estimator ──
    if W > 0:
        # Exonic RNA via strand subtraction
        exon_sense = spliced_sense[t] + unspliced_sense[t]
        exon_anti = spliced_anti[t] + unspliced_anti[t]
        exon_rna = max(0.0, (exon_sense - exon_anti) / denom)
        str_exon_density = exon_rna / max(1.0, L_exonic[t])

        # Intronic RNA via strand subtraction
        intron_rna = max(0.0, (intronic_sense[t] - intronic_anti[t]) / denom)
        if L_intronic[t] > 0:
            str_intron_density = intron_rna / L_intronic[t]
        else:
            str_intron_density = 0.0

        str_nrna = str_intron_density
        str_mrna = max(0.0, str_exon_density - str_intron_density)
    else:
        str_nrna = den_nrna
        str_mrna = den_mrna

    # ── Weighted combination ──
    final_nrna = W * str_nrna + (1.0 - W) * den_nrna
    final_mrna = W * str_mrna + (1.0 - W) * den_mrna

    # Evidence weight = total exonic + intronic fragments
    evidence = exon_total + intron_total

    return final_mrna, final_nrna, evidence
```

### 4.3 nrna_frac Prior from Hybrid Rates

```python
def compute_nrna_frac_from_rates(lambda_mrna, lambda_nrna):
    """Convert rate estimates to nascent fraction."""
    total = lambda_mrna + lambda_nrna
    if total > 0:
        return lambda_nrna / total
    else:
        return 0.5  # uninformative
```

### 4.4 nrna_frac Prior Cascade (refined)

The cascade uses the hybrid rate estimator at each level, aggregating
counts before computing rates (not averaging nrna_frac values):

```
Level 0 — Transcript:
    Aggregate: counts for transcript t only
    Compute: hybrid rates → nrna_frac_t
    Use if: evidence_t >= MIN_EVIDENCE (default 10)

Level 1 — TSS group:
    Aggregate: sum counts across all transcripts sharing TSS
    Compute: hybrid rates → nrna_frac_tss
    Use if: evidence_tss >= MIN_EVIDENCE

Level 2 — Locus-strand:
    Aggregate: sum counts across all same-strand transcripts in locus
    Compute: hybrid rates → nrna_frac_strand
    Use if: evidence_strand >= MIN_EVIDENCE

Level 3 — Fallback:
    nrna_frac = 0.5, κ = 2 (maximum entropy, very weak)
```

At each level, κ (confidence) equals the total fragment evidence at
that level:

```python
def cascade_nrna_frac_prior(t, tss_group, locus_strand_group, s, gdna_density):
    """Compute (alpha, beta) Beta prior for nrna_frac_t."""

    # Level 0: Transcript
    mrna_t, nrna_t, ev_t = compute_hybrid_rates(t, s, gdna_density)
    if ev_t >= MIN_EVIDENCE:
        nrna_frac = compute_nrna_frac_from_rates(mrna_t, nrna_t)
        kappa = ev_t
        return nrna_frac * kappa, (1 - nrna_frac) * kappa

    # Level 1: TSS group
    mrna_tss, nrna_tss, ev_tss = sum_hybrid_rates(tss_group, s, gdna_density)
    if ev_tss >= MIN_EVIDENCE:
        nrna_frac = compute_nrna_frac_from_rates(mrna_tss, nrna_tss)
        kappa = ev_tss
        return nrna_frac * kappa, (1 - nrna_frac) * kappa

    # Level 2: Locus-strand
    mrna_ls, nrna_ls, ev_ls = sum_hybrid_rates(locus_strand_group, s, gdna_density)
    if ev_ls >= MIN_EVIDENCE:
        nrna_frac = compute_nrna_frac_from_rates(mrna_ls, nrna_ls)
        kappa = ev_ls
        return nrna_frac * kappa, (1 - nrna_frac) * kappa

    # Level 3: Uninformative fallback
    return 1.0, 1.0  # Beta(1,1) = Uniform
```

### 4.5 Single-Exon Transcripts

Single-exon transcripts (L_intronic = 0) deserve special treatment:

- **Density estimator**: Cannot estimate D_intron for this transcript.
  We set D_intron = gdna_density (assume background only), making
  den_nrna = 0. This is conservative — nRNA from single-exon genes
  cannot be distinguished from mRNA or gDNA by density alone.

- **Strand estimator**: Still works normally — strand subtraction on
  exonic reads estimates total RNA regardless of intron presence.

- **nrna_frac prior**: Will cascade to TSS group or locus-strand level, where
  multi-exon transcripts provide intronic signal. If no neighboring
  multi-exon transcripts exist, nrna_frac falls to 0.5 (uninformative) and the
  EM resolves it from the data.

### 4.6 gDNA Initialization (θ_t and γ_ℓ)

The hybrid estimator also improves θ_t initialization:

```python
# θ_t init: gDNA-corrected total RNA
# Use the hybrid rates to estimate pure RNA contribution
mrna_t, nrna_t, _ = compute_hybrid_rates(t, s, gdna_density)
theta_t_init = (mrna_t + nrna_t) * L_exonic[t] + pseudocount

# γ_ℓ init: unchanged EB-shrunk gDNA density
gamma_l_init = locus_gdna_density[l] * locus_footprint[l]
```

## 5. Removing the Hard Strand Threshold

### 5.1 Current behavior

The pipeline currently enforces:
```python
if strand_specificity < 0.6:
    raise ValueError("Library is not strand-specific enough")
```

This rejects unstranded libraries entirely.

### 5.2 New behavior

With the unified hybrid model:

- **s ≥ 0.9**: W_strand ≥ 0.64. Strand signal dominates. Best accuracy.
- **0.7 ≤ s < 0.9**: W_strand ∈ [0.16, 0.64]. Both signals contribute.
  Slightly reduced accuracy but fully functional.
- **0.5 ≤ s < 0.7**: W_strand < 0.16. Density signal dominates.
  Functional but with reduced gDNA/nRNA separation for capture libraries.
- **s = 0.5**: W_strand = 0. Pure density mode. No strand information
  used. Still functional for standard total RNA-seq.

The hard threshold is removed. Instead, a warning is emitted for low
strand specificity:

```python
if s < 0.7:
    logger.warning(
        "Low strand specificity (%.2f). gDNA/nRNA estimates will "
        "rely primarily on coverage density. Results may be less "
        "accurate for capture-enriched libraries.", s
    )
```

### 5.3 Implications for fragment scoring

The strand model in `FragmentScorer` already uses the estimated strand
specificity to weight sense vs antisense likelihoods. For unstranded
libraries (s ≈ 0.5), the strand likelihood ratio approaches 1.0 for
all fragments, effectively disabling strand-based scoring. This is
correct behavior — no information means no penalty.

The fragment scoring changes required:
- Remove the hard `STRAND_DENOM_MIN` threshold that forces fragments
  to "strand-qualified" status
- Replace with a continuous strand weight that gracefully degrades
- Ensure that when s ≈ 0.5, strand scoring contributes negligible
  signal (likelihood ratio ≈ 1)

## 6. Data Flow

```
                    ┌──────────────────┐
                    │   BAM Scan       │
                    │   (C++ scanner)  │
                    └────────┬─────────┘
                             │
                    Per-transcript counts:
                    spliced_sense, spliced_anti,
                    unspliced_sense, unspliced_anti,
                    intronic_sense, intronic_anti
                             │
                    ┌────────▼─────────┐
                    │  Strand Model    │
                    │  (estimate s)    │
                    └────────┬─────────┘
                             │ s
                    ┌────────▼─────────┐
                    │  EB gDNA Prior   │
                    │  (intergenic     │
                    │   density)       │
                    └────────┬─────────┘
                             │ gdna_density, locus_gdna_density
                    ┌────────▼─────────┐
                    │  Hybrid Rate     │
                    │  Estimator       │
                    │  (per-transcript │
                    │   λ_mrna,λ_nrna) │
                    └────────┬─────────┘
                             │
                    ┌────────▼─────────┐
                    │  nrna_frac Prior Cascade │
                    │  t → TSS → LS →  │
                    │  fallback        │
                    └────────┬─────────┘
                             │ (α_t, β_t) per transcript
                    ┌────────▼─────────┐
                    │  θ_t init        │
                    │  (gDNA-corrected │
                    │   total RNA)     │
                    └────────┬─────────┘
                             │ θ_t, nrna_frac_t, γ_ℓ
                    ┌────────▼─────────┐
                    │  Locus EM        │
                    │  (linked model)  │
                    └──────────────────┘
```

## 7. Interaction with Linked RNA Model

This initialization plan feeds directly into the linked mRNA–nRNA EM
described in `linked_rna_model.md`:

- **θ_t** (total RNA) is initialized from gDNA-corrected hybrid rates
- **nrna_frac_t** (nascent fraction) is initialized from the nrna_frac prior cascade
  using hybrid density+strand estimation
- **γ_ℓ** (locus gDNA) is initialized from EB-shrunk intergenic density
  (unchanged)
- **(α_t, β_t)** Beta prior parameters for nrna_frac_t are passed to the C++
  EM solver for MAP estimation in the M-step

The EM then refines all parameters from the data, with the initialization
providing a biologically grounded starting point.

## 8. Implementation Phases

### Phase 1a: Hybrid rate estimator (Python)
- Implement `compute_hybrid_rates()` in `estimator.py`
- Requires: per-transcript strand-split counts (already in scan output),
  exonic/intronic lengths (already in `TranscriptGeometry`), gDNA density
  (already from EB prior), strand specificity s (already from strand model)
- Unit tests: verify rate isolation with known inputs

### Phase 1b: Revise nrna_frac prior cascade
- Replace current `compute_nrna_frac_priors()` (which uses raw unspliced/spliced
  ratios) with the hybrid-rate-based cascade
- Same 4-level hierarchy: transcript → TSS → locus-strand → fallback
- Unit tests: verify cascade behavior, single-exon handling

### Phase 1c: Revise θ_t initialization
- Use gDNA-corrected hybrid rates for θ_t init instead of raw unambig counts
- Ensures silent genes under gDNA don't get inflated initial abundances

### Phase 1d: Remove hard strand threshold
- Remove the `strand_specificity < 0.6` rejection in `pipeline.py`
- Add warning for s < 0.7
- Verify that fragment scoring degrades gracefully for low s

### Phase 2: Linked EM (C++ changes)
- As described in `linked_rna_model.md` Phase 2–3
- The initialization changes in Phase 1 are prerequisites

### Phase 3: Testing
- Simulated unstranded library → verify tool runs and produces reasonable
  output
- Simulated capture library with high exonic gDNA → verify strand
  estimator mitigates density bias
- Benchmark: compare initialization quality (nrna_frac_init vs EM-converged nrna_frac)
  across strand specificity levels

## 9. Required Count Arrays

The hybrid estimator needs per-transcript counts split by:
1. Splice status: spliced / unspliced-exonic / intronic
2. Strand: sense / antisense

The current scan phase already tracks these in `unambig_counts` with
columns indexed by `SpliceStrandSlot`:
- `SPLICED_ANNOT_SENSE`, `SPLICED_ANNOT_ANTI`
- `SPLICED_UNANNOT_SENSE`, `SPLICED_UNANNOT_ANTI`
- `UNSPLICED_SENSE`, `UNSPLICED_ANTI`

These provide spliced (annot + unannot) and unspliced exonic counts by
strand. For intronic counts by strand, the scan phase currently
accumulates `transcript_unspliced_sense` and
`transcript_unspliced_antisense` in the pre-EM loop (`scan.py`). These
are per-transcript intronic fragment counts (sense/antisense) that we
can use directly.

**No new scan-phase data collection is needed.** All required inputs
are already computed.

## 10. Configuration Parameters

| Parameter | Default | Location | Purpose |
|-----------|---------|----------|---------|
| `tss_window` | 200 | `EMConfig` | Fuzzy TSS grouping window (bp) |
| `nrna_frac_min_evidence` | 10 | `EMConfig` | Min fragments for hierarchy level |
| *(removed)* | ~~0.6~~ | ~~pipeline~~ | ~~Hard strand threshold~~ |

No new configuration parameters are introduced. The hybrid weighting
is fully determined by the estimated strand specificity s.



Implementation phases
Phase 1a: Hybrid rate estimator function ✅ (estimator.py: _compute_hybrid_nrna_frac_vec, compute_global_gdna_density)
Phase 1b: Revised nrna_frac prior cascade ✅ (estimator.py: compute_hybrid_nrna_frac_priors replaces compute_nrna_frac_priors)
Phase 1c: gDNA-corrected θ_t initialization (pending)
Phase 1d: Remove hard strand threshold ✅ (pipeline.py: warning at s<0.7 instead of rejection at s≤0.6)
Phase 2: Linked EM C++ changes (from linked_rna_model.md)
Phase 3: Unstranded + capture library testing
Ready to begin Phase 1a implementation on your signal.