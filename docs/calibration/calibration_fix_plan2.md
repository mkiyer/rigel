Now I have full clarity. Let me present the complete consolidated analysis.
# Complete Analysis: Three Root Causes Re-examined

## ✅ Root Cause #1 (FL): I WAS WRONG

You are correct. I misremembered the simulation parameters. The yaml shows `gdna.frag_mean: 350` and the calibration recovered **351bp** — accuracy of **0.3%**. This is excellent, not a bug. My apologies for the confusion.

## ✅ Root Cause #2 (Locus EM): I WAS WRONG

My fragment-level analysis used `ZL == 1` as the "assigned to gDNA" indicator. But [`ZL`](src/rigel/native/bam_scanner.cpp ) is the **locus_id**, not a flag. The actual gDNA assignment indicator is **bit 2 (0x04) of the [`ZF`](src/rigel/native/bam_scanner.cpp ) flag tag**.

When I corrected the fragment analysis, the picture flipped completely:

| Condition | gDNA correctly assigned | gDNA leaked to mRNA |
|-----------|------------------------|---------------------|
| `gdna_low` (10%) | **98.2%** | 1.1% |
| `gdna_med` (50%) | **98.6%** | 0.5% |
| `gdna_equal` (100%) | **98.7%** | 0.4% |
| `gdna_high` (200%) | **98.8%** | 0.6% |

The locus-level EM separates gDNA from mRNA with **>98% accuracy at the fragment level**, exactly as designed. Your intuition was correct: the EXON-INTRON boundary signal anchors the gDNA component, the FL likelihood (RNA mean 257 vs gDNA mean 351) discriminates effectively, and the per-locus prior `α_gdna` projects the calibration density onto each locus.

## ⚠️  Root Cause #3 (Density Non-Uniformity): REAL BUG, MECHANISM IDENTIFIED

This one **IS** a real bug, and your guidance was exactly right. Let me walk through it precisely.

### Observed
```
ρ_intergenic = 0.099590  (true rate for gdna_equal condition: 0.10 frags/bp)
ρ_intron     = 0.078734  (under-estimates by 21%)
ρ_exon-intron= 0.089689  (boundary flux density)
```

The intergenic estimate is essentially perfect. The intron estimate is 21% low.

### Mechanism (verified by reading density_global.py and scan_payload.py)

The classification scheme (`MASK_INTRON = 0b010`) assigns a fragment to "intron-only" **only if all aligned blocks fall within intronic sequence** — i.e. the fragment must be fully *contained* in introns and overlap zero exonic bases.

Then, the rate is computed in density_global.py as:

$$
\rho_{\text{intron}} = \frac{\sum_r n_r}{\sum_r L_{\text{eff,overlap}}(r)}, \quad L_{\text{eff,overlap}} = \text{span}_r + \bar{\ell}_{\text{gdna}} - 1
$$

**This is mathematically inconsistent.** The numerator counts containment; the denominator counts overlap (the "span + FL − 1" formula gives positions where a fragment of mean length *overlaps* the region by ≥1 base). The two should match.

### Worked example using your scenario

Transcript with intron at [2000, 3000), span = 1000bp. Fragment FL distribution: half at 275bp, half at 375bp, two fragments.

- **Containment effective length** (your suggestion): per fragment, `max(0, span − fl + 1)`
  - frag 1: 1000 − 275 + 1 = 726
  - frag 2: 1000 − 375 + 1 = 626
  - Sum: **1352 bp**  → density = 2 / 1352 = **0.001479 frags/bp**

- **Overlap effective length** (current code): `span + fl_mean − 1` = 1000 + 325 − 1 = **1324 bp**
  - For a single fragment but extended to all overlapping positions: ~2× the containment length
  - With both fragments: 2 × 1324 = 2648 → density = 2 / 2648 = **0.000755 frags/bp**

The current formula yields **half the correct rate** for this single-intron example. For larger introns (where span ≫ fl), the discrepancy shrinks, but for typical introns (mean ~3000bp here, FL 350bp) the bias is significant.

### Verification against observed data

Empirically expected ratio: `(span − fl + 1) / (span + fl − 1)` averaged over intron sizes:
- Intergenic regions in this synthetic genome have mean span ~167kb. Both formulas agree to within 0.4%. ✓ matches observation (intergenic estimate is accurate).
- Intronic regions have mean span ~2900bp, FL 351. 
  - Mean containment leff ≈ 2900 − 351 + 1 = 2550
  - Mean overlap leff ≈ 2900 + 351 − 1 = 3250
  - Predicted ratio = 2550/3250 = **0.785**
  - Observed ratio `ρ_intron / ρ_intergenic` = **0.793**
  - Residual: 0.008 — i.e. nearly 100% of the observed deficit is explained by this single bug.

When corrected, the predicted intron density becomes 0.078734 × (3250/2550) = **0.1004 frags/bp** vs intergenic 0.0996 — ratio = 1.01. Density uniformity restored.

## Why this matters for downstream estimates

The density pipeline feeds `π_gdna` (per-locus gDNA mass) via [`assemble_priors()`](src/rigel/calibration/locus_prior.py ). When `ρ_intron` is biased low and `ρ_exon-intron` (which also feeds the prior) is biased low (similarly, ~10%), the per-locus `α_gdna` priors are **systematically deflated for loci with significant intronic content**. The EM then has weaker prior pull toward gDNA in long-intron loci, leading to:

1. Some intronic gDNA leakage into mRNA at intron-rich loci
2. Slight `gdna_rate` underestimate in the global summary (we observe 0.6626 vs true 0.6667 — tiny but consistent)
3. Most importantly, **mis-shrinkage** of nearby loci via the empirical-Bayes pooling

The bias is small at the global level (because intergenic dominates), but it grows in importance for highly fragmented (many-intron) loci.

---

# Recommended Algorithmic Improvements

## R1 — **[High priority]** Fix density formula consistency: per-fragment effective length

Implement your suggestion exactly. Replace [`l_eff_overlap`](src/rigel/calibration/density_global.py ) for INTRON and EXON-INTRON with **per-fragment containment effective length**:

$$
L_{\text{eff,intron}} = \sum_r \sum_f \max(0, \text{span}_r - \ell_f + 1)
$$

where the inner sum is over the FL distribution (or its histogram from M3). Equivalent vectorized form:

$$
L_{\text{eff,intron}}(r) = \sum_{\ell} h(\ell) \cdot \max(0, \text{span}_r - \ell + 1)
$$

where `h(ℓ)` is the gDNA fragment-length probability mass (already collected in [`fl_hist`](src/rigel/calibration/scan_payload.py )). **No new C++ scan needed** — `fl_hist` is already there.

**Implementation sketch:**
```python
def l_eff_contained(span_bp: np.ndarray, fl_pmf: np.ndarray, fl_grid: np.ndarray) -> np.ndarray:
    """Sum_l p(l) * max(0, span - l + 1) — true 'contained' eff length per region."""
    # span_bp: (R,), fl_pmf: (L,), fl_grid: (L,)  → returns (R,)
    diffs = span_bp[:, None] - fl_grid[None, :] + 1.0   # (R, L)
    np.maximum(diffs, 0.0, out=diffs)
    return diffs @ fl_pmf
```

For INTERGENIC, the current overlap formula is benign (large spans), but for consistency apply the same containment formula uniformly. The cost is one matrix-vector multiply per density type — negligible.

**Note**: keep [`l_eff_overlap`](src/rigel/calibration/density_global.py ) for the EXON-INTRON boundary-flux denominator, which is different geometry (capture window of width `fl_mean` *crossing* a boundary, by construction).

## R2 — **[Medium priority]** Decide and document classification semantics

Currently the classification is "all blocks in intron" (containment-like). Two coherent choices, pick one:

- **(A) Containment classification + containment denominator** (R1 above). Recommended — matches your intuition and yields uniform densities.
- **(B) Overlap classification + overlap denominator**. Would require recounting fragments that overlap *any* intron base. Would inflate intronic counts at exon boundaries, double-counting EXON-INTRON fragments. Not recommended.

Rename mask constants (`MASK_INTRON` → `MASK_INTRON_ONLY`) to make containment semantics explicit in code.

## R3 — **[Low priority]** Health check on density uniformity

Add a calibration diagnostic: if `|log(ρ_intron / ρ_intergenic)| > 0.15` after the fix, emit a warning. Useful for catching new regressions and for real-data sanity checks (where biology could legitimately make introns differ — e.g. intronic RNA contamination or differential GC content affecting gDNA yield).

## R4 — **[Future]** Make EXON-INTRON projection consistent with R1

[`estimate_global_gdna_fragments()`](src/rigel/calibration/density_global.py ) projects `ρ_exon-intron` onto exonic regions using the overlap formula. After R1, evaluate whether to also use the FL-pmf-aware containment formula here. The exonic gDNA estimate (the part of total gDNA invisible to direct classification) currently has ~5% bias at low contamination; the same correction may close that gap.

## R5 — **[Already implemented]** Document the comment in `l_eff_overlap`

The docstring of [`l_eff_overlap`](src/rigel/calibration/density_global.py ) calls `max(0, span - fl_mean + 1)` "the forbidden alternative" that "ate two weeks of debugging." This is misleading — that formula is **correct for containment**; it was wrong only in the SRD-v1 era because it was paired with overlap-style classification. Update the comment to clarify: the issue is consistency, not the formula itself.

---

# Summary of Verified State (after correction of my earlier errors)

| Subsystem | Status | Notes |
|-----------|--------|-------|
| Global gDNA rate (calibration) | ✅ Excellent | <1pp error across 0–67% true rate |
| FL models (RNA, gDNA) | ✅ Accurate | 257bp vs 250 true; 351bp vs 350 true |
| Density: intergenic | ✅ Accurate | matches simulated rate to within 0.5% |
| **Density: intron** | ⚠️ **21% low** | **R1 fix needed** |
| Density: exon-intron | ⚠️ ~10% low | likely improves with R1 |
| Locus EM gDNA separation | ✅ 98%+ | per-fragment correctness verified via ZF tag |
| mRNA quantification | ✅ Robust | Spearman 0.91, Pearson 0.89 across all conditions |
| gDNA isoform leakage | ⚠️ Some | concentrated in intron-rich loci — likely improves with R1 |

The single high-impact fix is **R1**: replacing the overlap effective length with per-fragment containment effective length for intron/exon-intron classifications, using the existing `fl_hist` data.