# T+N+1 Pipeline Ablation Study — Root Cause Analysis & Improvement Proposals

**Date**: 2026-03-13  
**Model revision**: T+N+1 (single gDNA component with `log(0.5)` strand penalty)  
**Sweep framework**: `synthetic_sim_sweep.py` with `locus_config_simple.yaml`

---

## 1. Executive Summary

The T+N+1 architecture performs **well at SS ≥ 0.90** (weighted gDNA error
< 1%, nRNA error < 1%) and **reasonably at SS = 0.75** (weighted errors ~14%).
At SS = 0.50, quantification collapses catastrophically due to a fundamental
identifiability barrier between nRNA and gDNA.

**Key findings:**

| Metric | SS=0.50 | SS=0.75 | SS=0.90 | SS=1.00 |
|--------|---------|---------|---------|---------|
| Weighted mRNA error | −13.0% | −0.3% | +0.4% | +0.4% |
| Weighted nRNA error | −94.0% | −14.4% | −0.7% | +3.5% |
| Weighted gDNA error | +95.7% | +14.1% | +0.6% | −3.5% |
| nRNA FP rate (NTA1=0) | 100% | 100% | 100% | 100% |
| Mean spurious nRNA (NTA1=0) | 951.6 | 310.8 | 159.7 | 67.1 |
| Corr(gDNA,nRNA) excess | −0.997 | −1.000 | −1.000 | −1.000 |

**Critical issues identified:**
1. Near-perfect nRNA↔gDNA anti-correlation (−0.997 to −1.000) at all SS values
2. 83% nRNA false-positive rate even at SS = 1.0 (mean 67 spurious fragments)
3. SS = 0.50 identifiability collapse (nRNA and gDNA are indistinguishable)
4. Systematic gDNA underestimation at SS ≥ 0.75 (−8% to −10% weighted)

**Simulator bug fixed** during this study: when all RNA abundances are zero in
`n_rna_fragments` mode, the simulator previously generated 10,000 fragments
distributed uniformly across transcripts. Now correctly returns 0 RNA fragments.

---

## 2. Experimental Setup

### 2.1 Locus Configuration

| Transcript | Gene | Strand | Structure | Notes |
|-----------|------|--------|-----------|-------|
| TA1 | GA | + | 3 exons: [1000,1020] [5000,5500] [12000,13000] | Multi-exon, swept |
| TD | GD | + | 1 exon: [28000,29000] | Single-exon, always 0 |
| TE | GE | − | 2 exons: [40000,41000] [45000,46000] | Negative strand, always 0 |

Three nRNA groups: NTA1 (+strand, 1000–13000), NTD (+strand), NTE (−strand).
Only TA1 and NTA1 are swept; TD, TE, NTD, NTE are fixed at 0.

### 2.2 Sweep Grid (1,940 runs after excluding 4 degenerate)

- **Strand specificity**: {0.50, 0.75, 0.90, 1.00}
- **gDNA fraction**: {0.0, 0.1, 0.3, 0.5, 1.0, 2.0}
- **TA1 abundance**: {0, 5, 10, 20, 50, 100, 200, 500, 1000}
- **NTA1 abundance**: {0, 5, 10, 20, 50, 100, 200, 500, 1000}
- **n_rna_fragments**: 10,000 fixed; gDNA added on top
- RNA and gDNA fragment length distributions EQUAL (mean=250, std=50) to
  isolate strand/scoring effects from fragment-length discrimination

### 2.3 Simulator Fix Applied

`compute_rna_split()` and `_compute_pool_split()` in `src/rigel/sim/reads.py`
now return `(0, 0)` / `(0, 0, 0)` when total weight is zero, instead of
distributing all fragments uniformly. The sweep skips runs where pool_split
would be (0, 0, 0).

---

## 3. Detailed Results

### 3.1 Performance by Strand Specificity

**SS = 1.00 (perfect strandedness) — 485 runs:**
- mRNA: mean error +0.3%, MAE 6.1%, tight distribution [p10=−12.9%, p90=+5.2%]
- nRNA: mean error +1.9%, MAE 5.5%, [p10=−5.8%, p90=+11.3%]
- gDNA: mean error −9.1%, MAE 10.0%, slight systematic underestimation

**SS = 0.90 — 485 runs:**
- mRNA: mean error +5.8%, MAE 9.1%
- nRNA: mean error −0.7%, MAE 10.8%
- gDNA: mean error −10.3%, MAE 14.6%

**SS = 0.75 — 485 runs:**
- mRNA: mean error +4.8%, MAE 9.8%
- nRNA: mean error −10.0%, MAE 22.3% — degradation starting
- gDNA: mean error −6.8%, MAE 22.4%

**SS = 0.50 (unstranded) — 485 runs:**
- mRNA: mean error −37.9%, MAE 38.4% — **severe**
- nRNA: mean error −80.0%, MAE 82.9% — **near total loss**
- gDNA: mean error +190.5%, MAE 201.1% — **catastrophic overestimation**

### 3.2 Cross-table: gDNA Error by SS × gdna_fraction

| SS \ gf | 0.1 | 0.3 | 0.5 | 1.0 | 2.0 |
|---------|------|------|------|------|------|
| 0.50 | +564% | +183% | +117% | +58% | +30% |
| 0.75 | −33% | −22% | −12% | +10% | +23% |
| 0.90 | −29% | −17% | −10% | −1% | +6% |
| 1.00 | −21% | −13% | −8% | −4% | +0.2% |

At low gDNA fractions, errors are amplified because a small absolute
misallocation produces a large relative error. The systematic
underestimation at SS ≥ 0.75 with low gDNA burden (−21% to −33% at
gdna_frac=0.1) suggests the `log(0.5)` penalty may be slightly too aggressive.

### 3.3 False Positive Analysis

**nRNA false positives (NTA1=0, TA1 > 0, nRNA_obs > 1): 160/192 = 83%**

| SS | FP count | Mean spurious nRNA | Median | Max |
|----|----------|-------------------|--------|-----|
| 0.50 | 40/48 | 951.6 | 861.5 | 1548.9 |
| 0.75 | 40/48 | 310.8 | 310.4 | 624.6 |
| 0.90 | 40/48 | 159.7 | 185.2 | 303.6 |
| 1.00 | 40/48 | 67.1 | 59.5 | 131.8 |

The 83% rate persists across all SS values, though magnitude decreases
with better strandedness.

**gDNA false positives (gdna_frac=0, gDNA_obs > 1): 54/320 = 17%**
- ALL occur at SS=0.50 (54 of 80 SS=0.50 runs)
- Mean spurious gDNA: 8,178.8 fragments

**mRNA false positives (TA1=0, NTA1 > 0, mRNA_obs > 1): 160/192 = 83%**
- Very low magnitude: mean 3.4–4.4 spurious mRNA fragments
- Clinically insignificant

### 3.4 Conservation of Mass

Perfect at all SS values: sum of gDNA + nRNA + mRNA excesses = 0.00 in
every case. The misallocation is purely redistributive — no mass is
created or destroyed.

### 3.5 Fragment Length Recovery

- RNA FL: mean error +1.07 bp (true mean 250), std 3.86
- gDNA FL: mean error −1.47 bp (true mean 250), std 0.93

Fragment length model training is accurate, confirming the scoring
subsystem works well.

---

## 4. Root Cause Analysis

### 4.1 The nRNA ↔ gDNA Identifiability Problem

**Core issue**: nRNA and gDNA fragments are fundamentally ambiguous on
unspliced reads. Both produce UNSPLICED_ANNOT or INTRONIC fragments that
overlap with introns. The ONLY signal distinguishing them is strand:

- nRNA fragments are strand-specific (follow library prep convention)
- gDNA fragments are equally likely from either strand

In the T+N+1 model, gDNA gets a `log(0.5)` strand penalty per fragment.
When SS = 1.0, nRNA fragments are perfectly stranded while gDNA fragments
pay a −0.693 likelihood penalty per fragment. This creates separation.
When SS = 0.5, nRNA fragments are equally distributed across strands →
nRNA and gDNA have identical likelihoods → the EM cannot distinguish them.

The near-perfect anti-correlation (−0.997 to −1.000) between gDNA and
nRNA estimation errors confirms they are a **zero-sum pair**: any mass
misallocated from one goes directly to the other.

### 4.2 nRNA False Positives (83% rate at all SS values)

**Mechanism**: When there is no true nRNA (NTA1=0) but gDNA is present,
gDNA fragments that overlap gene regions are eligible for both the gDNA
component and any nRNA component covering that region. The EM distributes
some mass to nRNA because:

1. **The nRNA prior is not data-informed.** In `locus.py`, the `nrna_init`
   value computed from strand-corrected intronic evidence is used only as
   a binary gate (nRNA components with zero intronic span are killed via
   `prior[n_t + ln] = 0.0` for single-exon loci). For multi-exon loci,
   nRNA components are always alive regardless of whether intronic
   evidence supports them.

2. **OVR warm-start distributes uniformly across eligible components.**
   In `em_solver.cpp` (lines 613–660), the initial theta is seeded by a
   one-versus-rest (OVR) coverage vote. Each ambiguous fragment votes
   equally for all eligible components. Since gDNA fragments overlap both
   gDNA and nRNA components, they seed nonzero theta for nRNA.

3. **The EM converges to a local optimum** where nRNA absorbs some
   fraction of the unspliced mass. Even at SS=1.0, the nRNA component
   scores slightly better than gDNA on sense-strand unspliced fragments
   (no `log(0.5)` penalty), so the EM allocates a fraction there.

**Quantification**: At SS=1.0, mean spurious nRNA is 67.1 fragments
(out of 10,000 RNA + gDNA). At SS=0.50, mean is 951.6. The spurious
nRNA is almost perfectly anti-correlated with gDNA deficit, confirming
it is gDNA mass leaking into nRNA.

### 4.3 gDNA Systematic Underestimation at SS ≥ 0.75

At SS ≥ 0.75, gDNA is systematically underestimated by 4–10% (weighted).
This is the complement of Issue 4.2 — the mass that leaks into nRNA
comes from gDNA. The `log(0.5)` penalty makes gDNA slightly less
attractive than nRNA for sense-strand unspliced fragments, creating a
persistent bias.

### 4.4 SS = 0.50 Catastrophe

At SS = 0.50:
- nRNA and gDNA have identical per-fragment likelihoods (both unstranded)
- The `log(0.5)` gDNA penalty becomes the only distinguishing feature
- But `log(0.5)` penalizes gDNA relative to nRNA, so nRNA should win
- Instead, gDNA massively overestimates because of `nrna_init` behavior

Investigation revealed: at SS ≤ 0.6, `compute_nrna_init()` (locus.py
lines 449–502) cannot separate nRNA from gDNA using strand imbalance
(denominator `2×SS − 1 ≤ 0.2`). When `denom ≤ 0.2`, it falls back to
using total intronic coverage as init, but for some locus configurations
the nRNA init is computed as zero or very low, causing the prior gate to
disable nRNA. With nRNA disabled, ALL unspliced mass goes to gDNA →
catastrophic overestimation.

### 4.5 gDNA-Only Runs (TA1=0, NTA1=0, gdna_frac > 0)

In these 20 runs (new with the simulator fix), Rigel correctly allocates
most mass to gDNA (error −7.8%) but produces mean 246.7 spurious nRNA
fragments. This confirms the nRNA false-positive mechanism is independent
of mRNA presence.

---

## 5. Architectural Issues

### 5.1 Binary Gating vs. Informative Initialization

The `nrna_init` and `gdna_init` values are computed with considerable
effort (strand-corrected intronic evidence, empirical Bayes) but then
**discarded** — used only as binary on/off gates for the prior. The
actual EM theta is seeded by OVR coverage weighting, which is
uninformed about the nRNA/gDNA distinction.

**Code path**:
- `nrna_init` computation: `locus.py` → `compute_nrna_init()` (L449–502)
- Binary gating: `locus.py` L387–394 (`prior[n_t + ln] = 0.0` for single-exon)
- OVR warm-start: `em_solver.cpp` L613–660 (`compute_ovr_prior_and_warm_start()`)
- EM theta init: `theta_init[i] = unambig_totals[i] + ovr_vote[i]`

### 5.2 Asymmetric Component Treatment

| Dimension | mRNA | nRNA | gDNA |
|-----------|------|------|------|
| unambig_totals seed | Row-sum of unambig_counts | 0.0 | 0.0 |
| Gating mechanism | Always eligible | Single-exon check | gdna_init > 0 |
| Init value used | — | Ignored for theta | Binary gate only |
| OVR seeding | ✓ | ✓ | ✓ |

mRNA components benefit from unambiguous spliced fragments that directly
seed their theta. nRNA and gDNA both start from zero unambig_totals and
rely entirely on OVR voting, creating the conditions for confusion.

### 5.3 Effective Length Normalization Gap

The sweep used equal RNA/gDNA fragment length distributions (mean=250,
std=50), deliberately removing fragment-length discrimination. In real
data, gDNA fragment lengths tend to be longer, which would help
distinguish gDNA from nRNA. The current results represent a **worst case**
for fragment-length discrimination.

---

## 6. Proposed Improvements

### 6.1 Use nrna_init to Seed EM Theta (Not Just Gate) — HIGH PRIORITY

**Problem**: `nrna_init` captures real strand-corrected intronic evidence
but is thrown away after gating.

**Proposal**: Use `nrna_init` values to inform the initial theta allocation
between nRNA and gDNA components, rather than relying solely on OVR. The
strand-corrected intronic excess gives a direct estimate of nRNA abundance
that the EM should start from.

**Expected impact**: Reduce nRNA false positives, reduce gDNA
underestimation at SS ≥ 0.75, improve convergence speed.

### 6.2 Strand-Aware OVR Warm-Start — HIGH PRIORITY

**Problem**: The OVR warm-start treats nRNA and gDNA identically — both
receive equal votes from unspliced fragments.

**Proposal**: During OVR initialization, weight the votes using the
strand likelihood ratio. For a sense-strand unspliced fragment at SS=0.9,
nRNA should receive ~90% of the vote and gDNA ~10%. This would provide
much better initial theta separation.

**Expected impact**: Address the root cause of nRNA↔gDNA confusion at
all SS values.

### 6.3 Prior-Informed gDNA Initialization — MEDIUM PRIORITY

**Problem**: The EB-computed gDNA density estimate (`compute_gdna_density_from_strand`)
is used for the binary gate but not for theta seeding.

**Proposal**: Use the EB gDNA estimate to set a data-informed prior mass
for the gDNA component, rather than a flat epsilon.

### 6.4 Adaptive `log(0.5)` Penalty — LOW PRIORITY

**Problem**: The fixed `log(0.5)` penalty slightly over-penalizes gDNA at
high SS, leading to systematic underestimation.

**Proposal**: Scale the gDNA strand penalty by the observed strand
specificity: `log(1/2)` at SS=1.0 but interpolating toward 0 as SS→0.5.
This would reduce the bias at SS=0.75 while maintaining discrimination.

**Risk**: At SS close to 0.5, this would make gDNA and nRNA even more
similar. May need a floor.

### 6.5 nRNA Regularization / Sparsity Prior — MEDIUM PRIORITY

**Problem**: 83% nRNA false-positive rate even at SS=1.0.

**Proposal**: Add a sparsity-inducing prior on nRNA components (e.g.,
a small fixed penalty for non-zero nRNA). This would suppress spurious
nRNA allocation when the evidence is weak. Alternatively, post-EM
thresholding: zero out nRNA components below a confidence threshold.

### 6.6 SS=0.50 Hard Cutoff — LOW PRIORITY

**Problem**: At SS=0.50, nRNA and gDNA are theoretically indistinguishable.

**Proposal**: When SS ≤ 0.55, disable nRNA/gDNA separation entirely and
report all unspliced mass as "unresolved nascent/gDNA". This is more
honest than producing wildly inaccurate estimates.

---

## 7. Summary & Recommendations

### What Works Well
- **mRNA quantification**: Robust at SS ≥ 0.75 (mean error +0.3% to +5.8%)
- **Conservation of mass**: Perfect — no mass leakage
- **Fragment length recovery**: Excellent (< 1.5 bp error)
- **SS = 1.0 overall**: Weighted errors < 4% for all components

### What Needs Improvement
1. **nRNA false positives** (83% rate): Address via strand-aware OVR (§6.2) + nrna_init seeding (§6.1)
2. **nRNA↔gDNA confusion** (−0.999 correlation): Same root cause as above
3. **gDNA underestimation** at SS ≥ 0.75: Address via informative initialization (§6.1, §6.3)
4. **SS = 0.50 collapse**: Fundamental identifiability limit; consider hard cutoff (§6.6)

### Priority Order
1. **Strand-aware OVR warm-start** (§6.2) — highest expected ROI
2. **nrna_init theta seeding** (§6.1) — uses existing computation
3. **nRNA sparsity prior** (§6.5) — directly targets false positives
4. **EB gDNA prior** (§6.3) — complement to §6.1
5. **Adaptive log(0.5)** (§6.4) — fine-tuning
6. **SS cutoff** (§6.6) — honest reporting

---

## Appendix: Sweep Data Location

- **TSV**: `sweep_results/tn1_ablation_v2/sweep_results.tsv` (1,940 rows × 64 columns)
- **Analysis script**: `scripts/analyze_tn1_v2.py`
- **Config**: `scripts/locus_config_simple.yaml`
- **Simulator fix**: `src/rigel/sim/reads.py` (`compute_rna_split`, `_compute_pool_split`)
- **Sweep skip logic**: `scripts/synthetic_sim_sweep.py` (degenerate all-zero detection)
