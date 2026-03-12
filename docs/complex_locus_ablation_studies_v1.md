# Complex Locus Ablation Studies v1

## Objective

Identify and diagnose accuracy failures in Rigel's EM initialization and prior
framework by running controlled synthetic experiments on a complex locus with
overlapping isoforms, convergent genes, single-exon genes, and varying levels
of nascent RNA (nRNA) and genomic DNA (gDNA) contamination.

The ultimate goal is to **simplify** what we believe is an over-engineered,
over-parameterized prior framework (27 EMConfig fields, 5 distinct prior
mechanisms) while improving accuracy.

---

## Experimental Design

### Complex Locus

A synthetic locus containing 10 transcripts spanning diverse structural
challenges:

| Transcript | Strand | Structure | Purpose |
|------------|--------|-----------|---------|
| TA1 | + | 4-exon, long | Overlapping isoform family (Gene A) |
| TA2 | + | 3-exon | Subset of TA1 exons |
| TA3 | + | 2-exon | Minimal overlap with TA1 |
| TA4 | + | 2-exon, short | Shortest A-family isoform |
| TB1 | - | 3-exon | Antisense overlap with Gene A |
| TB2 | - | 2-exon | Antisense overlap with Gene A |
| TCF | + | 1-exon | Single-exon, convergent pair |
| TCR | - | 3-exon | Multi-exon, convergent with TCF |
| TD  | + | 1-exon | Single-exon, isolated |
| TE  | - | 2-exon | Negative control (unexpressed) |

8 nRNA groups: NTA1, NTA2, NTA34, NTB12, NTCR, NTE (user-defined) +
nrna_TCF, nrna_TD (auto-derived for single-exon genes).

### Controlled Variables

- **Fragment length**: Equal for RNA and gDNA (mean=250, std=50) — isolates
  non-FL issues
- **RNA fragments**: Fixed at 10,000 — gDNA added *on top* so RNA truth
  values don't change with contamination level
- **Read length**: 150 bp, paired-end

### Sweep Dimensions

| Dimension | Values | Count |
|-----------|--------|-------|
| Expression pattern | 16 correlated patterns (see below) | 16 |
| Strand specificity (SS) | 0.5, 0.75, 0.9, 1.0 | 4 |
| gDNA fraction | 0.0, 0.1, 0.3, 0.5 | 4 |
| TB1 expression | 0, 500 | 2 |
| TB2 expression | 0, 500 | 2 |

**Total: 16 × 4 × 4 × 2 × 2 = 1,024 runs**

### 16 Expression Patterns

| # | Name | Category | Key Feature |
|---|------|----------|-------------|
| 0 | TA1_only_no_nrna | Pure mRNA | Single isoform, no nRNA |
| 1 | TA4_only_no_nrna | Pure mRNA | Short isoform only |
| 2 | mixed_A_no_nrna | Pure mRNA | All 4 A isoforms |
| 3 | convergent_CF_CR_no_nrna | Pure mRNA | Convergent gene pair |
| 4 | TA1_mod_nrna | Moderate nRNA | 200 frags/nRNA |
| 5 | TA4_mod_nrna | Moderate nRNA | Short isoform + nRNA |
| 6 | mixed_A_mod_nrna | Moderate nRNA | All A + moderate nRNA |
| 7 | convergent_mod_nrna | Moderate nRNA | Convergent + nRNA |
| 8 | TA1_heavy_nrna | Heavy nRNA | 1000 frags/nRNA |
| 9 | mixed_A_heavy_nrna | Heavy nRNA | All A + heavy nRNA |
| 10 | TCF_heavy_nrna | Heavy nRNA | Single-exon + heavy nRNA |
| 11 | TB_heavy_nrna | Heavy nRNA | Antisense gene heavy nRNA |
| 12 | everything_moderate | Complex | All components active |
| 13 | TA1+TCF+heavy_nrna | Complex | Mixed heavy nRNA |
| 14 | TA3+TCR+heavy_nrna | Complex | Cross-strand heavy nRNA |
| 15 | nrna_only_no_mrna | Complex | **Stress test**: only nRNA |

---

## Results Summary

All 1,024 runs completed. Analysis covers 10 sections. Key findings below.

### 1. Total mRNA Accuracy

| nRNA Level | SS=0.5 | SS=0.75 | SS=0.9 | SS=1.0 | Notes |
|------------|--------|---------|--------|--------|-------|
| None, gdna=0 | 2.5% | **0.0%** | **0.0%** | **0.0%** | Near-perfect baseline |
| None, gdna=0.5 | 9.7% | 5.4% | 3.4% | 1.5% | Modest gDNA leakage |
| Moderate, gdna=0 | 36.8% | 10.2% | **12.9%** | 1.1% | nRNA causes mRNA errors even without gDNA |
| Moderate, gdna=0.5 | 45.8% | 21.9% | 16.4% | 16.1% | Severe |
| Heavy, gdna=0 | 45.0% | 5.8% | 5.5% | 5.2% | Better than moderate (!) |
| Heavy, gdna=0.5 | 64.5% | 64.0% | 38.4% | 23.1% | gDNA interaction severe |
| Complex, gdna=0 | 48.6% | 52.7% | 47.4% | 7.4% | Catastrophic at SS<1.0 |
| Complex, gdna=0.5 | 67.4% | 67.6% | **106.8%** | 102.5% | Max single run: **864%** |

The paradox that heavy nRNA (1000 frags) recovers *better* than moderate nRNA
(200 frags) at SS≥0.75 without gDNA is a critical clue — see Root Cause 2.

### 2. nRNA Accuracy

| nRNA Level | SS=0.5 | SS=0.75 | SS=0.9 | SS=1.0 |
|------------|--------|---------|--------|--------|
| Moderate, gdna=0 | **-83%** | -18% | -14% | -0.1% |
| Moderate, gdna=0.5 | -90% | -14% | -19% | -24% |
| Heavy, gdna=0 | -75% | **-0.4%** | -0.4% | 0.0% |
| Heavy, gdna=0.5 | -71% | +6% | +6% | +3% |
| Complex, gdna=0 | -88% | -20% | -13% | -0.3% |
| Complex, gdna=0.5 | -94% | -20% | -18% | -23% |

Massive nRNA underestimation at moderate levels (-83% at SS=0.5). Heavy nRNA
recovers well at SS≥0.75 without gDNA, but gDNA causes systematic nRNA
overestimation (+2-6%).

### 3. gDNA Accuracy

| nRNA Level | SS=0.5 | SS=0.75 | SS=0.9 | SS=1.0 | Pattern |
|------------|--------|---------|--------|--------|---------|
| None | +13-26% | -11 to -31% | -7 to -22% | -5 to -14% | Reasonable |
| Moderate | +139-685% | **+0.4 to +60%** | +9 to +66% | +15 to +84% | nRNA→gDNA misattribution |
| Heavy | +136-719% | -20 to -38% | -16 to -33% | -9 to -29% | Underestimate |
| Complex | +174-815% | +17 to +162% | +12 to +111% | +19 to +127% | Severe overestimate |

At SS=0.5, gDNA is massively overestimated across all nRNA levels — nRNA mass
is being dumped into gDNA. At SS≥0.75 with moderate nRNA, gDNA is
*overestimated* while nRNA is *underestimated* — a clear misattribution.

### 4. Per-Transcript Errors

- **TA1**: Systematic -7% to -20% underestimate (mass leaks to TA4)
- **TA4**: +17% to +28% overestimate at SS≥0.75 (absorbs mass from longer
  isoforms)
- **TB1**: -34% at SS=0.5, near-perfect at SS≥0.75
- **TCR**: **+45% to +77%** overestimate at SS≥0.75 with gDNA — convergent gene
  confusion

### 5. False Positives (Unexpressed Transcripts)

| Transcript | Total FP Cases | Worst Mean FP (frags) | Key Trigger |
|------------|---------------|----------------------|-------------|
| TD | 587 | 147.8 | gDNA + any nRNA |
| TE | 637 | 51.4 | gDNA (any pattern) |
| TCR | 447 | **544.0** | nRNA patterns at SS≥0.75 |
| TCF | 404 | 105.5 | gDNA only |
| TB1 | 368 | 38.0 | nRNA + gDNA |
| TA4 | 325 | 30.4 | When TA1 only expressed |
| TA2 | 385 | 7.5 | When TA1/TA4 expressed |
| TA1 | 294 | 69.8 | When TA4 expressed + gDNA |

Single-exon genes (TD, TCF) and the convergent gene TCR are the most
vulnerable to false positives.

### 6. Worst Cases

All top-20 worst runs are **pattern 15** (nRNA-only, no mRNA expressed):

| SS | gDNA frac | mRNA error | nRNA error | gDNA error |
|----|-----------|-----------|-----------|-----------|
| 0.9 | 0.5 | **+864%** | -1262 | -258 |
| 0.9 | 0.3 | +763% | -1715 | +373 |
| 1.0 | 0.5 | +756% | -1378 | +47 |
| 0.75 | 0.5 | +375% | -1245 | +247 |

The EM cannot push the nRNA fraction η to 1.0 (all nRNA, no mRNA), so exonic
nRNA reads leak into mRNA as false positives.

### 7. nRNA Stress Test (Pattern 15)

At SS=0.5 with TB genes expressed: nRNA drops to **0.1 fragments** while gDNA
inflates to **9,990+**. The nRNA component is being killed entirely.

At SS≥0.75 without gDNA and no TB: near-perfect (nRNA=9,947, mRNA=53). But
adding TB1=500 at SS=0.75: mRNA balloons to 1,179 (expected 266), nRNA
collapses to 7,283 (expected 9,734).

### 8. Convergent Gene Accuracy

**Pattern 3** (TCF+TCR convergent pair, no nRNA):
- SS≥0.75, gdna=0: **Near-perfect** (< 1% error for both TCF and TCR)
- SS=0.5, gdna=0.5: TCF -55%, TCR -44% (strand confusion at low SS)
- SS≥0.75, gdna=0.5: Slight overestimate (+2-7%) from gDNA leakage

**Pattern 7** (convergent + moderate nRNA):
- TCF at SS=0.5: **-100%** (completely missed)
- TCR at SS≥0.75 with lower expression: **+140% to +179%** (nRNA reads
  misattributed)
- TCR at SS≥0.75 with higher expression: +30-35% (less severe)

### 9. gDNA Contamination Effect on mRNA

| gDNA frac | SS=0.5 | SS=0.75 | SS=0.9 | SS=1.0 |
|-----------|--------|---------|--------|--------|
| 0.0 | -812 | +249 | +244 | +9 |
| 0.1 | -994 | +525 | +507 | +428 |
| 0.3 | -1116 | +646 | +548 | +453 |
| 0.5 | -1264 | +696 | +554 | +475 |

At SS≥0.75, gDNA contamination drives **mRNA overestimation** (+430 to +700
fragment excess), compensated by nRNA underestimation. The mass balance shows
nRNA absorbs ~2× the compensation needed.

### 10. Antisense Overlap Effect

Mean Gene A absolute error with TB off vs on:

| gDNA frac | SS=0.75 TB off | SS=0.75 TB on | SS=0.9 TB off | SS=0.9 TB on |
|-----------|---------------|--------------|--------------|-------------|
| 0.0 | 129 | 103 | 129 | 103 |
| 0.1 | 140 | 121 | 136 | 118 |
| 0.3 | 140 | 117 | 133 | 112 |
| 0.5 | 123 | 101 | 123 | 97 |

Surprisingly, antisense expression *reduces* Gene A error at SS≥0.75 —
possibly because it provides better strand-based disambiguation. At SS=0.5,
the pattern reverses (TB on increases error).

---

## Root Cause Analysis

### RC1: nRNA Components Killed at SS ≤ 0.6 (Catastrophic)

**Symptom**: Pattern 15 at SS=0.5 with antisense TB genes: nRNA gets 0.1
fragments while gDNA inflates to 9,990+.

**Code path**: `compute_nrna_init()` has a hard threshold `2·SS − 1 ≤ 0.2` →
returns all zeros. Then `build_locus_em_data()` sees `nrna_init == 0` → zeros
the nRNA prior → component is permanently dead. All nRNA reads must go
somewhere → the strand symmetry penalty (which rewards balanced strand
distributions) funnels everything into gDNA at SS=0.5.

**Impact**: Unstranded or poorly-stranded libraries (~50% of public data) lose
**all** nRNA resolution. This is a hard gate with no recovery path.

### RC2: 3-Tier nRNA EB Over-Shrinks Moderate nRNA Signals (Severe)

**Symptom**: Moderate nRNA (200 frags/nRNA): underestimated -14% to -83%,
gDNA overestimated +60-85%. Heavy nRNA (1000 frags) recovers much better.

**Code path**: `compute_nrna_frac_priors()` performs hierarchical shrinkage:

$$\hat{\eta}_{nrna} = \frac{ev_{nrna} \cdot \eta_{nrna} + \kappa_{locus} \cdot \hat{\eta}_{ls}}{ev_{nrna} + \kappa_{locus}}, \quad \hat{\eta}_{ls} = \frac{ev_{ls} \cdot \eta_{ls} + \kappa_{global} \cdot \bar{\eta}}{ev_{ls} + \kappa_{global}}$$

When per-nRNA evidence is small (moderate nRNA), shrinkage toward the
locus-strand and global means dominates. The MoM κ estimates are noisy with
few observations. The 10 tuning parameters create an opaque parameter surface.

**Key insight**: Heavy nRNA recovers better because it has enough signal to
overcome the EB shrinkage bias. This is the EB working *against* the data.

### RC3: nRNA Exonic Mass Leaks to mRNA — η Cannot Reach 1 (Catastrophic)

**Symptom**: Pattern 15 (pure nRNA, no mRNA): up to +864% mRNA false positive
error. This is the single worst failure mode.

**Code path**: The hierarchical M-step computes:

$$\eta_{MAP} = \frac{c_{nrna} + \alpha - 1}{c_{nrna} + c_{mrna} + \alpha + \beta - 2}$$

Even when $c_{mrna}$ should be zero, the Beta prior prevents η from reaching
1.0. The OVR warm start seeds mRNA with nonzero mass via coverage-weighted
shares. Once mRNA has mass, the EM gets stuck.

**Fundamental**: Exonic nRNA reads have identical signature to mRNA reads —
an inherent identifiability problem that the prior framework should mitigate,
not amplify.

### RC4: Convergent Gene nRNA → TCR Overestimation (Severe)

**Symptom**: Pattern 7 (convergent + moderate nRNA): TCR gets +140% to +179%
error at SS≥0.75 with lower expression levels.

**Code path**: neg-strand nRNA intronic reads look like TCR mRNA to the EM.
The nRNA deduplication by unique genomic span doesn't prevent cross-strand
nRNA mass from being misattributed. At lower expression levels, the
signal-to-noise is worse and the shrinkage prior prevents proper nRNA
attribution.

### RC5: Single-Exon Gene False Positives from gDNA (Moderate)

**Symptom**: TD: 587 FPs (highest count), TE: 637 FPs, TCF: 404 FPs (only
with gDNA > 0).

**Code path**: gDNA reads uniformly cover single-exon genes. The OVR warm
start distributes gDNA mass to all eligible components. Single-exon genes
lack splice junctions to disambiguate.

**Partially inherent**: Single-exon genes are genuinely ambiguous with gDNA.
But the OVR prior compounds the problem by seeding these components with
nonzero initial mass.

### RC6: gDNA → mRNA Leakage at SS ≥ 0.75 (Moderate)

**Symptom**: gDNA fragments cause mRNA overestimation (+430 to +700 Δ at
SS≥0.75), compensated by nRNA underestimation (-680 to -1110 Δ).

**Code path**: The gDNA EB hierarchy underestimates gDNA's contribution
(3-tier shrinkage may over-regularize), surplus mass flows into mRNA/nRNA.
The nRNA component absorbs the compensation because it's the most flexible
component.

---

## Assessment of the 5 Prior Mechanisms

| Mechanism | Params | Assessment |
|-----------|--------|------------|
| Flat Dirichlet α = 0.01 | 1 | **OK** — small pseudocount, minimal harm |
| OVR coverage-weighted γ = 1.0 | 1 | **Harmful** — seeds FPs by distributing gDNA mass to unexpressed genes |
| nRNA 3-tier EB (MoM κ) | 10 | **Over-engineered** — shrinkage hurts moderate nRNA; 10 params for marginal benefit |
| gDNA 3-tier EB | 4 | **Over-engineered** — reference-level aggregation adds complexity; global intergenic density would suffice |
| Component eligibility gating | 0 (implicit) | **Too aggressive** — kills nRNA at SS ≤ 0.6, making it irreversible |
| Strand symmetry penalty | 2 | **Context-dependent** — helpful at high SS, destructive at low SS |

---

## Phased Implementation Plan

The plan is ordered by **expected accuracy impact × simplicity of change**.
Each phase includes:
1. A detailed design document (before implementation)
2. Implementation
3. Re-running the 1,024-run complex locus ablation to measure improvement

### Phase 1: Fix nRNA Component Gating (RC1)

**Priority**: Highest — this is a catastrophic bug that kills nRNA at SS ≤ 0.6.

**Difficulty**: Low — localized change in `compute_nrna_init()` and
`build_locus_em_data()`.

**What to change**:
- Remove the hard `2·SS − 1 ≤ 0.2` threshold in `compute_nrna_init()` that
  zeros all nRNA initialization
- At low SS, fall back to total intronic coverage (both strands) as the nRNA
  init signal, rather than returning zeros
- In `build_locus_em_data()`, do not permanently zero nRNA priors when
  `nrna_init == 0` — instead, use a small but nonzero floor so the EM can
  discover nRNA signal from the data
- Keep the single-exon nRNA zeroing (genuinely no intronic span to measure)

**Expected impact**: Restore nRNA resolution for unstranded/poorly-stranded
libraries. Should dramatically improve SS=0.5 accuracy across all sections.

**Validation**: Re-run 1,024 ablation. Expect SS=0.5 nRNA accuracy to improve
from -83% to something reasonable. gDNA overestimation at SS=0.5 should drop.

---

### Phase 2: Eliminate 3-Tier nRNA Empirical Bayes (RC2)

**Priority**: Very high — the hierarchical shrinkage hurts moderate nRNA and
creates 10 opaque parameters.

**Difficulty**: Moderate — touches `compute_nrna_frac_priors()` (priors.py)
and the hierarchical M-step (em_solver.cpp).

**What to change**:
- Replace the 3-tier shrinkage pipeline (global → locus-strand → per-nRNA)
  with a direct per-nRNA estimate
- Use the existing `_compute_hybrid_nrna_frac_vec()` density+strand estimator
  at the per-nRNA level (this is already computed — currently it's the input
  to the EB, now it becomes the final estimate)
- Replace the MoM-estimated `κ_global` and `κ_locus` shrinkage with a simple
  flat Beta prior — either Beta(1,1) (uniform) or Beta(0.5, 0.5) (Jeffreys)
- Keep `kappa_nrna` as the single tuning parameter for prior strength, or
  eliminate it entirely with a flat prior

**Parameters removed**: 9 of 10 nRNA-related params (`kappa_global`,
`kappa_locus`, `kappa_min`, `kappa_max`, `kappa_fallback`, `kappa_min_obs`,
`mom_min_evidence_global`, `mom_min_evidence_locus`, `kappa_nrna`) → 1 or 0

**Expected impact**: Moderate nRNA accuracy should improve substantially
(currently -14 to -83%). May resolve the paradox where heavy nRNA recovers
better than moderate. Reduces parameter surface from 10 to 1.

**Validation**: Re-run 1,024 ablation. Compare moderate vs heavy nRNA
accuracy. Expect the gap to narrow.

---

### Phase 3: Replace OVR Prior with Flat + gDNA-Aware Warm Start (RC3, RC5)

**Priority**: High — OVR seeds false positives, especially for unexpressed
transcripts and single-exon genes.

**Difficulty**: Moderate — touches `compute_ovr_prior_and_warm_start()` in
em_solver.cpp and the prior construction in `build_locus_em_data()`.

**What to change**:
- Remove the OVR coverage-weighted prior term
  (`γ · coverage_i / n_ambig`) — use flat α only
- Redesign warm start: initialize mRNA proportional to unambiguous counts
  (splice-informed signal only), not total coverage
- For gDNA warm start: use the per-locus strand-based estimate directly
- For nRNA warm start: use the per-nRNA density-based estimate directly
- This prevents gDNA mass from seeding unexpressed genes at initialization

**Parameters removed**: 1 (`prior_gamma`)

**Expected impact**: Should reduce false positives across TD (587→fewer), TE
(637→fewer), TCF (404→fewer). Should reduce the +864% mRNA error in
pattern 15 by not seeding mRNA with gDNA-derived coverage.

**Validation**: Re-run 1,024 ablation. Focus on false positive counts and
pattern 15 mRNA error.

---

### Phase 4: Scale Strand Symmetry Penalty by SS (RC1, RC6)

**Priority**: Medium — the penalty is helpful at high SS but destructive at
low SS.

**Difficulty**: Low — localized change in the M-step of em_solver.cpp.

**What to change**:
- Scale the strand symmetry κ by strand specificity:
  $\kappa_{eff} = \kappa \cdot (2 \cdot SS - 1)^2$
- At SS=0.5: penalty is zero (no strand information available)
- At SS=1.0: full penalty strength
- Fold the `strand_symmetry_pseudo` regularization into the κ scaling,
  eliminating one parameter

**Parameters removed**: 1 (`strand_symmetry_pseudo`)

**Expected impact**: At SS=0.5, prevents the penalty from funneling RNA mass
into gDNA. At SS=0.75, reduces over-penalization. Works synergistically with
Phase 1 (nRNA gating fix).

**Validation**: Re-run 1,024 ablation. Focus on SS=0.5 and SS=0.75 gDNA
accuracy.

---

### Phase 5: Simplify gDNA Initialization to Per-Locus Estimate (RC6)

**Priority**: Medium — the 3-tier gDNA EB adds complexity for marginal
benefit.

**Difficulty**: Moderate — touches `compute_eb_gdna_priors()` (locus.py) and
fan-out/aggregation logic.

**What to change**:
- Replace the global → reference → locus 3-tier EB hierarchy with direct
  per-locus estimation
- Per-locus gDNA density: use the strand-asymmetry formula
  $G = 2(anti \cdot SS - sense \cdot (1-SS))/(2SS-1)$ computed from local
  unspliced counts
- Floor at global intergenic density (already available from
  `compute_global_gdna_density()`)
- Remove MoM κ estimation for reference and locus levels
- `gdna_init = max(strand_estimate, intergenic_floor) × exonic_bp`

**Parameters removed**: 4 (`gdna_kappa_ref`, `gdna_kappa_locus`,
`gdna_mom_min_evidence_ref`, `gdna_mom_min_evidence_locus`)

**Expected impact**: Cleaner gDNA estimates. Removes reference-level
aggregation that doesn't help in practice (in our sweep, all loci share one
reference anyway). Reduces parameter count.

**Validation**: Re-run 1,024 ablation. Compare gDNA accuracy columns across
all nRNA levels.

---

### Phase 6: Address nRNA-mRNA Identifiability (RC3)

**Priority**: Medium-low — this is the hardest problem and partially inherent.

**Difficulty**: High — requires careful changes to the hierarchical M-step in
C++ and potentially the component model itself.

**What to change** (requires detailed investigation):
- Allow η to reach 1.0 (or very close) when evidence strongly supports pure
  nRNA — the current Beta(α,β) MAP prevents this
- Consider Jeffreys prior Beta(0.5, 0.5) which is less restrictive than the
  current informative prior
- Investigate whether splice-junction information can explicitly distinguish
  mRNA from nRNA (spliced reads → mRNA, unspliced exonic reads → either)
- Consider a simpler redistribution: allocate exonic mass to mRNA only in
  proportion to spliced evidence, with remainder going to nRNA

**Expected impact**: Should reduce the catastrophic +864% mRNA error in
pattern 15. But this is fundamentally a hard identifiability problem — exonic
nRNA reads are indistinguishable from mRNA reads without splice junctions.

**Validation**: Re-run 1,024 ablation. Focus on pattern 15 and moderate nRNA
mRNA accuracy.

---

### Phase 7: Convergent Gene nRNA Handling (RC4)

**Priority**: Low — affects a specific structural configuration.

**Difficulty**: High — may require changes to nRNA span deduplication and
cross-strand attribution logic.

**What to change** (requires detailed investigation):
- When nRNAs from different strands overlap, intronic reads from one strand's
  nRNA can be misattributed to the other strand's mRNA
- May need strand-aware nRNA component indexing
- May need to modify how nRNA intronic reads on one strand are excluded from
  mRNA components on the opposite strand

**Expected impact**: Should reduce TCR overestimation (+140-179%) in
convergent gene scenarios. Niche but not negligible.

**Validation**: Re-run 1,024 ablation. Focus on Section 8 (convergent gene
accuracy).

---

### Summary of Expected Parameter Reduction

| Phase | Params Removed | Running Total |
|-------|---------------|---------------|
| Start | — | 27 |
| Phase 1 (nRNA gating) | 0 (removes magic number) | 27 |
| Phase 2 (nRNA EB) | 9 | 18 |
| Phase 3 (OVR prior) | 1 | 17 |
| Phase 4 (strand penalty SS) | 1 | 16 |
| Phase 5 (gDNA EB) | 4 | 12 |
| Phase 6 (nRNA-mRNA ident.) | 0 | 12 |
| Phase 7 (convergent genes) | 0 | 12 |
| **Final** | **15** | **12** |

---

## Appendix: Files and Outputs

- **Config**: `scripts/complex_locus_config.yaml`
- **Analysis script**: `scripts/analyze_complex_locus.py`
- **Sweep results**: `/Users/mkiyer/Downloads/rigel_runs/complex_locus/sweep_results.tsv`
- **Full analysis output**: `/Users/mkiyer/Downloads/rigel_runs/complex_locus/analysis.txt`
- **Sweep log**: `/Users/mkiyer/Downloads/rigel_runs/complex_locus/sweep.log`

---

## Appendix: EMConfig Parameters (Current 27)

| Parameter | Default | Phase |
|-----------|---------|-------|
| `prior_alpha` | 0.01 | Keep |
| `prior_gamma` | 1.0 | Phase 3: Remove |
| `mode` | "map" | Keep |
| `iterations` | 1000 | Keep |
| `convergence_delta` | 1e-6 | Keep |
| `prune_threshold` | 0.1 | Keep |
| `confidence_threshold` | 0.95 | Keep |
| `nrna_frac_kappa_global` | None (auto) | Phase 2: Remove |
| `nrna_frac_kappa_locus` | None (auto) | Phase 2: Remove |
| `nrna_frac_kappa_nrna` | 5.0 | Phase 2: Remove or simplify |
| `nrna_frac_mom_min_evidence_global` | 50.0 | Phase 2: Remove |
| `nrna_frac_mom_min_evidence_locus` | 20.0 | Phase 2: Remove |
| `nrna_frac_kappa_min` | 2.0 | Phase 2: Remove |
| `nrna_frac_kappa_max` | 200.0 | Phase 2: Remove |
| `nrna_frac_kappa_fallback` | 5.0 | Phase 2: Remove |
| `nrna_frac_kappa_min_obs` | 20 | Phase 2: Remove |
| `gdna_kappa_ref` | None (auto) | Phase 5: Remove |
| `gdna_kappa_locus` | None (auto) | Phase 5: Remove |
| `gdna_mom_min_evidence_ref` | 50.0 | Phase 5: Remove |
| `gdna_mom_min_evidence_locus` | 30.0 | Phase 5: Remove |
| `strand_symmetry_kappa` | 6.0 | Keep (Phase 4: scale by SS) |
| `strand_symmetry_pseudo` | 50.0 | Phase 4: Remove (fold into κ) |
| `frag_len_*` | various | Keep (unrelated) |
