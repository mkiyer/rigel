# Calibration Algorithm V2: Data-Driven Length-Weighted EM

## Motivation

The V1 calibration algorithm suffers from structural problems:

1. **Brittle fallback.** A hard `if n_eligible < 100` threshold switches from a 50-iteration EM with three convergent signals to a zero-iteration algebraic rule. This creates a discontinuous jump in logic and cannot be patched into robustness.

2. **Equal-weight Gaussian M-step.** A 5bp boundary region with 1 read exerts the same influence on $\mu_G$ and $\sigma^2_G$ as a 50,000bp intronic region. Under a Poisson model, the 5bp region's density estimate has variance $\sim 1/n$, making it noise-dominated.

3. **Arbitrary thresholds.** `GDNA_INIT_MIN_REGIONS = 100`, `GDNA_INIT_DENSITY_PERCENTILE = 0.10`, `SS > 0.55` — all hard gates that create discontinuities. The underlying statistical signals (strand LLR, density, FL) degrade continuously; the algorithm should too.

4. **Multimapper contamination.** Region counts include multimapped fragments (fractionally weighted), which pollute density estimates in regions with heavy multimapping. Calibration should operate on unique-mapper evidence only for reliability.

5. **Single-exon transcript conflation.** Unspliced-only regions from single-exon transcripts are conflated with intergenic gDNA regions during seeding. These regions have genuine RNA expression but no splice evidence, causing them to incorrectly boost the gDNA seed.

The V2 algorithm addresses all of these by: (1) length-weighting the M-step, (2) deriving priors empirically from the data itself, (3) using strand decomposition for seeding with continuous SS handling, and (4) enforcing strict region eligibility.

---

## Algorithm Overview

```
1. Compute summary statistics (per-region counts, density, strand balance)
2. Build eligibility mask (exclude multimappers, antisense-overlap, zero-count)
3. Classify region roles (intergenic, exonic, intronic, single-exon-only)
4. Strand decomposition → per-region gDNA/RNA estimates (where applicable)
5. Seed partition (intergenic → gDNA, spliced → RNA, strand-decomposed → informed)
6. Compute empirical Bayes priors from seed regions
7. EM iteration with length-weighted MAP M-step
8. Output per-region posteriors γ, global density λ_G, κ, FL model
```

No fallback. No hard thresholds. The EM always runs.

---

## Phase 1: Length-Weighted Gaussian M-Step

### Problem

The current M-step at `_m_step()` ([calibration.py:895–952](src/rigel/calibration.py#L895)) fits density Gaussians using region-level weights:

$$\mu_G = \frac{\sum_r \gamma_r \cdot x_r}{\sum_r \gamma_r}, \quad \sigma^2_G = \frac{\sum_r \gamma_r \cdot (x_r - \mu_G)^2}{\sum_r \gamma_r}$$

where $x_r = \log(n_r / L_r + \epsilon)$. Each region is one observation of log-density, regardless of its physical size. A 5bp region with 1 read (density $= 0.2\text{ frags/bp}$) gets the same vote as a 50,000bp region with 5 reads (density $= 0.0001\text{ frags/bp}$).

### Fix

Weight by $\gamma_r \cdot L_r$:

$$\mu_G = \frac{\sum_r \gamma_r \cdot L_r \cdot x_r}{\sum_r \gamma_r \cdot L_r}, \quad \sigma^2_G = \frac{\sum_r \gamma_r \cdot L_r \cdot (x_r - \mu_G)^2}{\sum_r \gamma_r \cdot L_r}$$

Same for the RNA component with $(1 - \gamma_r)$.

### Justification

Under a Poisson model, $n_r \sim \text{Poi}(\lambda L_r)$. The Fisher information about $\lambda$ from region $r$ is proportional to $L_r$. Longer regions provide more precise density estimates and should exert proportionally more influence on the density Gaussians.

### Impact

Large intergenic/intronic regions (50+ kbp) anchor $\mu_G$ and $\sigma^2_G$. Tiny boundary regions (5–50bp) from the boundary-sweep algorithm have negligible influence. This stabilizes the Gaussians even with very few regions.

### Implementation

Modify `_m_step()` and the initial Gaussian computation in `calibrate_gdna()` (lines 1156–1167):

```python
# Current (region-level)
w_g_sum = float(w_g.sum())
mu_g = float(np.sum(w_g * ld) / w_g_sum)

# Proposed (length-weighted)
rl = region_length[eligible]
wL_g = w_g * rl
wL_g_sum = float(wL_g.sum())
mu_g = float(np.sum(wL_g * ld) / wL_g_sum)
var_g = max(float(np.sum(wL_g * (ld - mu_g) ** 2) / wL_g_sum), _EPS)
```

The `if w_g_sum >= 1.0` safety net changes to `if wL_g_sum > _EPS` (any non-zero weighted length). This safety net is further superseded by Phase 4 (Bayesian priors).

**Files:** `src/rigel/calibration.py` — `_m_step()`, initial Gaussian computation in `calibrate_gdna()`

---

## Phase 2: Strict Region Eligibility

### Problem

The current eligibility mask is minimal: `(n_total > 0) & (region_length > 0)`. This admits:

1. **Multimapper-contaminated regions.** Region counts include all fragments regardless of NH ([resolve_context.h:1395](src/rigel/native/resolve_context.h#L1395)). Regions in multimapper hotspots (paralogs, repeats) have inflated counts from fractionally-allocated multimappers. These inflate density estimates and corrupt the gDNA/RNA separation.

2. **Antisense-overlap regions.** Regions with `tx_pos AND tx_neg` (both-strand transcription) have `gene_strand == 0`. The strand decomposition cannot operate here (sense vs antisense is undefined when RNA comes from both strands). The strand LLR already returns 0 for these regions, but they still contribute to density estimation where they act as outliers.

3. **Single-exon-only exonic regions.** A region that is exonic but has zero splice signal is either: (a) a single-exon transcript's exon, or (b) an exon with insufficient coverage to observe splicing. In either case, the region has genuine RNA signal but looks like "unspliced-only" to the seeding algorithm. Currently, these regions can incorrectly seed the gDNA component.

### Fix

Define three region classification tiers based on annotation and data context. These are not hard exclusions from the EM — they control which roles a region can play during seeding and strand decomposition.

**Region flags** (computed once from `stats` and `region_df`):

| Flag | Condition | Meaning |
|------|-----------|---------|
| `is_intergenic` | `NOT tx_pos AND NOT tx_neg` | No annotated transcripts. Reads here are likely gDNA or unannotated expression. |
| `is_exonic` | `exon_pos OR exon_neg` | Overlaps annotated exon. Has expected RNA expression. |
| `has_splice` | `n_spliced > 0` | Definitive RNA evidence. Hard constraint: $\gamma = 0$. |
| `is_single_strand` | `gene_strand != 0` | Unambiguous strand context. Strand decomposition is applicable. |
| `is_antisense_overlap` | `tx_pos AND tx_neg` | Both-strand transcription. Strand decomposition cannot operate. |

**Eligibility for calibration EM** (the `eligible` mask):

```python
eligible = (
    (stats["n_total"] > 0)
    & (stats["region_length"] > 0)
)
```

This remains broad — all regions with data participate in the EM posterior computation. But the *seeding* and *strand decomposition* use narrower subsets (see Phase 3).

**Multimapper exclusion.** This requires a code change in the C++ `RegionAccumulator`. Currently, `accumulate()` adds ALL fragments to `counts[]` regardless of `is_unique`. The fix: accumulate unique-mapper counts into a separate 4-column array (`counts_unique[]`), and pass this to Python as `region_counts_unique`. Calibration uses `region_counts_unique`; the full `region_counts` (with multimappers) is retained for downstream locus_gamma aggregation if needed.

Alternatively, track a per-region multimapper fraction and use it post-hoc. But the cleanest approach is a separate unique-only counter in C++.

**Implementation:**

1. **C++ change** ([resolve_context.h](src/rigel/native/resolve_context.h)): Add `counts_unique` array alongside `counts`. Only accumulate into `counts_unique` when `is_unique == true`. Expose both arrays to Python via the BamScanner result dictionary.

2. **Python change** ([pipeline.py](src/rigel/pipeline.py)): Extract `region_counts_unique` from the scan result and pass to `calibrate_gdna()` instead of (or alongside) `region_counts`.

3. **Calibration change** ([calibration.py](src/rigel/calibration.py)): `compute_region_stats()` operates on the unique-only counts.

**Files:** `src/rigel/native/resolve_context.h`, `src/rigel/native/bam_scanner.cpp`, `src/rigel/pipeline.py`, `src/rigel/calibration.py`

---

## Phase 3: Strand-Decomposed Seeding

### Problem

The current seeding approach compares unspliced region densities against spliced region densities using an empirical CDF quantile ([calibration.py:686–705](src/rigel/calibration.py#L686)). This has two issues:

1. **Requires an arbitrary density percentile** (`GDNA_INIT_DENSITY_PERCENTILE = 0.10`).
2. **Requires a minimum seed size** (`GDNA_INIT_MIN_REGIONS = 100`) — directly causing the fallback.
3. **Single-exon transcripts corrupt the gDNA seed.** Expressed single-exon genes produce unspliced-only reads at high density, which are then compared against the spliced-region eCDF. If their density is above the 10th percentile (which it usually is, since they're expressed), they avoid the gDNA seed. But if their density happens to be low, they incorrectly enter the gDNA seed.

### Fix: Analytical Strand Decomposition

For each region with single-strand annotation context (`gene_strand != 0`) and $SS > 0.5$, analytically decompose the unspliced reads into gDNA and RNA components using the known strand specificity.

**The model.** Under R1-antisense convention (dUTP/TruSeq Stranded), for a region on the sense strand of a gene:

- RNA sense reads (R1 on antisense strand) appear as neg-strand alignments
- RNA antisense reads (R1 on sense strand) appear as pos-strand alignments
- gDNA reads are symmetric (equal pos and neg)

Let $k$ = sense-strand reads (gene-relative), $n$ = total unspliced reads. The strand decomposition variables already exist in the codebase (see `compute_sense_fraction` at [calibration.py:280–312](src/rigel/calibration.py#L280)):

- For `gene_strand == +1`: $k_{\text{sense}} = n_{\text{unspliced}} - n_{\text{pos}}$
- For `gene_strand == -1`: $k_{\text{sense}} = n_{\text{pos}}$

The system of equations:

$$k_{\text{sense}} = n_{\text{RNA}} \cdot SS + n_{\text{gDNA}} \cdot 0.5$$
$$n - k_{\text{sense}} = n_{\text{RNA}} \cdot (1 - SS) + n_{\text{gDNA}} \cdot 0.5$$

Solving for $n_{\text{gDNA}}$:

$$\hat{n}_{\text{gDNA}} = \frac{(n - k_{\text{sense}}) - n \cdot (1 - SS)}{SS - 0.5} = \frac{n \cdot SS - k_{\text{sense}}}{SS - 0.5}$$

$$\hat{n}_{\text{RNA}} = n - \hat{n}_{\text{gDNA}}$$

Clamp both to $[0, n]$.

### Behavior Across the SS Spectrum

The denominator $(SS - 0.5)$ controls the decomposition's precision:

| SS | Denominator | Behavior |
|----|-------------|----------|
| 1.0 | 0.5 | Perfect decomposition. Antisense reads = gDNA/2 exactly. |
| 0.9 | 0.4 | Good decomposition. Small amplification of noise (×2.5). |
| 0.7 | 0.2 | Moderate decomposition. Noise amplified ×5. |
| 0.55 | 0.05 | Poor decomposition. Noise amplified ×20. Nearly useless. |
| 0.50 | 0.0 | Degenerate. Division by zero. Decomposition undefined. |

The decomposition is not a binary on/off — it **continuously degrades** as $SS \to 0.5$. The variance of $\hat{n}_{\text{gDNA}}$ under Poisson sampling is:

$$\text{Var}[\hat{n}_{\text{gDNA}}] = \frac{n \cdot SS \cdot (1-SS)}{(SS - 0.5)^2}$$

At $SS = 0.9$, $n = 100$: $\text{Var} \approx 56$, so $\text{SD} \approx 7.5$ — decent.
At $SS = 0.55$, $n = 100$: $\text{Var} \approx 9900$, so $\text{SD} \approx 99.5$ — garbage.

**Rather than applying a hard SS threshold**, we use the decomposition's own uncertainty to weight its influence. The seeding algorithm should trust strand-decomposed estimates less as $SS \to 0.5$.

### Seeding Algorithm

**Inputs:** Region stats with unique-mapper counts (Phase 2), strand specificity $SS$, region metadata.

**Step 1: Classify regions.**

```python
has_splice    = n_spliced > 0
is_intergenic = ~tx_pos & ~tx_neg     # No annotated transcripts
is_single_strand = gene_strand != 0   # Strand decomposition applicable
is_antisense_overlap = tx_pos & tx_neg  # Both-strand transcription
unspliced_only = ~has_splice & (n_total > 0) & eligible
```

**Step 2: Hard constraints.**

- `has_splice` → $\gamma = 0$ (expressed, hard constraint). These regions have definitive RNA evidence (spliced reads). This constraint is never relaxed.

**Step 3: Intergenic seed.**

- `is_intergenic & (n_total > 0)` → $\gamma = 1$ (gDNA seed). Intergenic regions have no annotated transcripts. Any reads here are overwhelmingly gDNA (or unannotated expression at very low levels). These are the most reliable gDNA seed regions.

**Step 4: Strand decomposition for annotated single-strand regions.**

For regions with `is_single_strand & unspliced_only & ~is_intergenic`:

Compute $\hat{n}_{\text{gDNA},r}$ and $\hat{n}_{\text{RNA},r}$ via the analytical decomposition.

Compute per-region gDNA fraction: $\hat{f}_{G,r} = \hat{n}_{\text{gDNA},r} / n_r$.

Confidence-weight by the decomposition's SNR. Define:

$$w_r = \frac{(SS - 0.5)^2 \cdot n_r}{SS \cdot (1 - SS)}$$

This is the inverse-variance weight (proportional to $1/\text{Var}[\hat{f}_G]$). Regions with high $w_r$ have reliable decompositions; low $w_r$ regions are uncertain.

Use the weighted decomposition for seeding:

- If $\hat{f}_{G,r} > 0.8$ and $w_r > 4$ (roughly 2σ significance): gDNA seed ($\gamma = 1$)
- If $\hat{f}_{G,r} < 0.2$ and $w_r > 4$: RNA seed ($\gamma = 0$)
- Otherwise: ambiguous ($\gamma = 0.5$, let the EM decide)

The "$w_r > 4$" threshold is a signal-to-noise gate, not an arbitrary region count. It means "the decomposition is at least 2σ significant." With $SS = 0.9$ and $n = 20$: $w = 0.16 \cdot 20 / 0.09 = 35.6$ → significant. With $SS = 0.55$ and $n = 20$: $w = 0.0025 \cdot 20 / 0.2475 = 0.20$ → not significant. The decomposition's own uncertainty naturally gates itself — no hard SS threshold needed.

**Step 5: Regions without strand decomposition.**

For regions with `gene_strand == 0` (antisense-overlap or intergenic with data) that are `unspliced_only`: these cannot be strand-decomposed. Set $\gamma = 0.5$ (ambiguous). The EM's density and FL signals will classify them.

**Step 6: Degenerate case — no spliced regions at all.**

If `has_splice.sum() == 0` (no spliced reads anywhere — pure gDNA or extreme low-expression):

- Intergenic regions still seed gDNA ($\gamma = 1$).
- All other regions get $\gamma = 0.5$.
- The Bayesian priors (Phase 4) stabilize the M-step.

### Impact

- **No `GDNA_INIT_DENSITY_PERCENTILE`.** The density eCDF comparison is replaced by strand decomposition + intergenic seeding.
- **No `GDNA_INIT_MIN_REGIONS`.** The minimum seed size guarantee is unnecessary because Bayesian priors (Phase 4) stabilize the M-step even with zero gDNA seed regions.
- **Single-exon transcripts are handled correctly.** Exonic regions on single-exon genes have `tx_pos OR tx_neg = True` (not intergenic). If they have `gene_strand != 0` and their strand decomposition shows high RNA fraction, they seed as RNA. They never incorrectly enter the gDNA seed.
- **Antisense-overlap regions are excluded from strand decomposition.** Regions with `gene_strand == 0` due to both-strand transcription get $\gamma = 0.5$ and are classified by density/FL only.
- **Multimappers excluded** (Phase 2). The decomposition operates on clean unique-mapper counts.

**Files:** `src/rigel/calibration.py` — rewrite `_seed_initial_partition()`, add `_strand_decompose()` helper, remove `density_percentile` and `min_gdna_regions` parameters

---

## Phase 4: Empirical Bayes MAP Priors

### Problem

The current M-step crashes when one component has fewer than 1.0 effective regions of weight (`if w_g_sum >= 1.0`). The fallback (`mu_g, var_g = mu_overall, var_overall`) is a band-aid. With length-weighting (Phase 1), this safety net changes units but doesn't solve the fundamental problem: the MLE Gaussian cannot be fit to zero data.

### Fix: Data-Derived Normal-Inverse-Gamma Priors

Instead of hardcoded universal constants, derive the prior means and variances from the seed regions computed in Phase 3.

**Prior mean computation (after seeding):**

$$\mu_{0,G} = \frac{\sum_{r \in \text{gDNA seed}} L_r \cdot x_r}{\sum_{r \in \text{gDNA seed}} L_r}$$

$$\mu_{0,R} = \frac{\sum_{r \in \text{expressed seed}} L_r \cdot x_r}{\sum_{r \in \text{expressed seed}} L_r}$$

These are the length-weighted mean log-densities of the gDNA and RNA seed regions, respectively. They are empirical — derived entirely from the data, not hardcoded. The gDNA prior mean is anchored by intergenic regions (which are plentiful even in small genomes). The RNA prior mean is anchored by spliced regions.

**Prior variance:**

$$\sigma^2_{0,G} = \max\left(\frac{\sum_{r \in \text{gDNA seed}} L_r \cdot (x_r - \mu_{0,G})^2}{\sum_{r \in \text{gDNA seed}} L_r}, \; 1.0\right)$$

The floor of 1.0 prevents the variance prior from over-constraining (log-density variance across regions is typically 2–10 in real data).

**Prior strength $\kappa_0$:**

$$\kappa_0 = \alpha_{\text{prior}} \cdot \sum_{r \in \text{eligible}} L_r$$

where $\alpha_{\text{prior}}$ is a small fraction (default 0.01 = 1% of total eligible basepairs). This is the only new hyperparameter.

Scaling behavior:
- **Large dataset** (human genome, ~100k regions, ~3Gbp eligible): $\kappa_0 \approx 30\text{M bp}$. M-step sees ~3Gbp·γ of weighted data. Prior is $\sim 1\%$ — a whisper.
- **Small dataset** (test scenario, 5 regions, ~12kbp eligible): $\kappa_0 \approx 120\text{ bp}$. M-step sees ~6kbp·γ. Prior is $\sim 2\%$ — enough to prevent collapse, too weak to fight the data.

The prior weight scales with the data, keeping the prior-to-data ratio constant regardless of genome size.

**MAP M-step update:**

$$\mu_G^{\text{MAP}} = \frac{\kappa_0 \cdot \mu_{0,G} + \sum_r \gamma_r L_r x_r}{\kappa_0 + \sum_r \gamma_r L_r}$$

$$\sigma^{2,\text{MAP}}_G = \frac{\kappa_{0,\sigma} \cdot \sigma^2_{0,G} + \sum_r \gamma_r L_r (x_r - \mu_G^{\text{MAP}})^2}{\kappa_{0,\sigma} + \sum_r \gamma_r L_r}$$

where $\kappa_{0,\sigma}$ is the variance prior strength. We can set this equal to $\kappa_0$ (same scaling) or use a small constant (e.g., $\kappa_{0,\sigma} = 2 \cdot \bar{L}$ where $\bar{L}$ is the mean region length — equivalent to "2 virtual regions worth of variance evidence"). The choice matters little when data is abundant.

**Degenerate case — no gDNA seed regions:**

If the intergenic + strand-decomposed seeding produces zero gDNA seed regions (e.g., pure RNA sample from a small targeted panel with no intergenic sequence):

- $\mu_{0,G}$ cannot be computed from seed. Fall back to: $\mu_{0,G} = \mu_{0,R} - 3.0$ (gDNA density is assumed ~20× lower than RNA density). This is a weak, conservative assumption.
- The prior strength $\kappa_0$ ensures the M-step does not crash. The EM will converge with $\gamma \approx 0$ everywhere if there is no gDNA signal.

Similarly, if no expressed seed (no spliced regions): $\mu_{0,R} = \mu_{0,G} + 3.0$.

### Impact

- The `if w_g_sum >= 1.0` / `mu_overall, var_overall` safety net is deleted. The prior prevents degenerate M-step updates.
- The entire algebraic fallback block is deleted (Phase 5).
- The EM M-step always produces well-defined $\mu$ and $\sigma^2$, regardless of the number of regions.

**Files:** `src/rigel/calibration.py` — modify `_m_step()` to accept prior parameters, compute priors during initialization. `src/rigel/config.py` — add `prior_weight_fraction: float = 0.01` to `CalibrationConfig`.

---

## Phase 5: Delete Fallback and Hard Thresholds

### Deletions

With Phases 1–4 in place, the EM is stable for any number of regions. The following code is deleted:

1. **The algebraic fallback block** in `calibrate_gdna()` ([calibration.py:1063–1133](src/rigel/calibration.py#L1063)) — ~70 lines. The `if n_eligible < min_gdna_regions:` branch and all its contents.

2. **`GDNA_INIT_MIN_REGIONS`** constant ([calibration.py:70](src/rigel/calibration.py#L70)).

3. **`GDNA_INIT_DENSITY_PERCENTILE`** constant ([calibration.py:69](src/rigel/calibration.py#L69)).

4. **`min_gdna_regions`** and **`density_percentile`** from `CalibrationConfig` ([config.py:174–175](src/rigel/config.py#L174)) and from `calibrate_gdna()` signature.

5. **The `SS > 0.55` gate** — this was inside the deleted fallback block. The strand LLR natively degrades to zero as $SS \to 0.5$ (Beta-Binomial/Binomial LLR is identically zero at $SS = 0.5$), and the strand decomposition's confidence-weighting naturally disables it at low SS. No hard threshold needed anywhere.

6. **The `if w_g_sum >= 1.0` / `mu_overall, var_overall` fallback** in `_m_step()` — replaced by MAP prior regularization.

### Result

A single code path. `calibrate_gdna()` always runs: stats → seed → EM loop → output. Five regions or five hundred thousand — identical logic.

**Files:** `src/rigel/calibration.py`, `src/rigel/config.py`, `src/rigel/pipeline.py` (remove parameter forwarding)

---

## Detailed Data Flow

```
                         ┌──────────────────────┐
                         │  BAM Scan (C++)       │
                         │  ┌─────────────────┐  │
                         │  │ RegionAccumulator│  │
                         │  │ counts_unique[]  │──── unique-mapper 4-col counts
                         │  │ fl_table         │──── unique, unspliced, single-region FL
                         │  └─────────────────┘  │
                         └──────────────────────┘
                                    │
                    ┌───────────────┴───────────────┐
                    ▼                               ▼
          region_counts_unique              fl_table (unchanged)
                    │
                    ▼
        ┌───────────────────────┐
        │ compute_region_stats()│ ← operates on unique counts only
        │ (n_pos, n_neg,       │
        │  n_spliced, n_total, │
        │  density, gene_strand│
        │  tx_pos, tx_neg,     │
        │  exon_pos, exon_neg) │
        └───────────┬──────────┘
                    │
             ┌──────┴──────┐
             ▼             ▼
      Region flags    Strand decomposition
      (intergenic,    (n_gDNA, n_RNA per region
       single_strand,  where gene_strand ≠ 0)
       has_splice,
       antisense_ovl)
             │             │
             └──────┬──────┘
                    ▼
        ┌───────────────────────┐
        │ Seed partition         │
        │  • has_splice → γ=0   │
        │  • intergenic → γ=1   │
        │  • strand_decomp →    │
        │    confidence-weighted │
        │  • ambiguous → γ=0.5  │
        └───────────┬──────────┘
                    │
                    ▼
        ┌───────────────────────┐
        │ Empirical Bayes priors │
        │  μ₀_G from gDNA seeds │
        │  μ₀_R from RNA seeds  │
        │  κ₀ = 0.01 × Σ L_r   │
        └───────────┬──────────┘
                    │
                    ▼
        ┌───────────────────────┐
        │ EM Loop (MAP)          │
        │  E-step: density +    │
        │    strand + FL LLRs   │
        │  M-step: length-wtd   │
        │    MAP Gaussians      │
        │  κ update, FL rebuild │
        │  convergence: |Δπ|    │
        └───────────┬──────────┘
                    │
                    ▼
        ┌───────────────────────┐
        │ GDNACalibration output │
        │  region_posteriors γ   │
        │  gdna_density_global   │
        │  kappa_strand          │
        │  gdna_fl_model         │
        │  mixing_proportion     │
        └───────────────────────┘
```

---

## Region Classification Reference

The following table summarizes how different region types are handled:

| Region Type | `tx_pos` | `tx_neg` | `exon_pos/neg` | `gene_strand` | `n_spliced` | Treatment |
|-------------|----------|----------|----------------|---------------|-------------|-----------|
| Intergenic | F | F | F | 0 | 0 | gDNA seed ($\gamma = 1$) |
| Intronic (sense-only) | T | F | F | +1 | 0 | Strand decomposition; seed by confidence |
| Exonic (sense-only, multi-exon gene) | T | F | T | +1 | >0 | Expressed seed ($\gamma = 0$, hard) |
| Exonic (sense-only, single-exon gene) | T | F | T | +1 | 0 | Strand decomposition; seed by confidence. Will NOT incorrectly seed gDNA. |
| Exonic (antisense-overlap) | T | T | T | 0 | varies | $\gamma = 0.5$, EM classifies via density/FL. No strand decomposition. |
| Intronic (antisense-overlap) | T | T | F | 0 | 0 | $\gamma = 0.5$, EM classifies via density/FL |
| Zero-count | any | any | any | any | 0 | Not eligible. $\gamma = \pi$ (prior). |

---

## Handling of Edge Cases

### Pure RNA (zero gDNA)

Intergenic regions have zero or near-zero reads. Strand decomposition shows $\hat{n}_{\text{gDNA}} \approx 0$ in annotated regions. gDNA seed is sparse (intergenic regions with very few reads). Prior $\mu_{0,G}$ is anchored by these sparse intergenic regions (correct — they represent the true gDNA background level near zero). EM converges to $\gamma \approx 0$ everywhere. Result: $\lambda_G \approx 0$, $\pi \approx 0$.

### Pure gDNA (zero RNA)

No spliced regions exist. Expressed seed is empty. $\mu_{0,R}$ falls back to $\mu_{0,G} + 3.0$. Intergenic regions seed gDNA. Strand decomposition shows $\hat{n}_{\text{RNA}} \approx 0$ everywhere (reads are strand-symmetric). EM converges to $\gamma \approx 1$ everywhere. Result: $\lambda_G = $ global density, $\pi \approx 1$.

### Tiny genome (5–10 regions)

Bayesian priors stabilize the M-step. Intergenic regions (at least some exist in any genome) seed gDNA. Spliced regions (if any) seed RNA. Length-weighting ensures that even with 5 regions, the largest ones anchor the Gaussians. The EM converges to a reasonable (if imprecise) classification.

### Unstranded library (SS = 0.5)

Strand decomposition produces $\hat{n}_{\text{gDNA}} = 0/0$ (undefined). Confidence weight $w_r = 0$. These regions get $\gamma = 0.5$ (ambiguous). Strand LLR in the E-step is identically zero. Classification relies entirely on density and FL signals. Intergenic regions still seed gDNA (no strand needed). The algorithm degrades gracefully — no special case, no `if` statement.

### Weakly stranded library (SS = 0.55)

Strand decomposition produces noisy estimates (variance amplification ×20). Confidence weight $w_r$ is very low for most regions (need $n > 100$ to reach significance). Most regions get $\gamma = 0.5$. Intergenic regions still seed gDNA reliably. The EM relies primarily on density + FL, supplemented by weak strand signal. This is correct — a barely-stranded library should provide barely-any strand information.

### Single-exon transcripts

Exonic regions on single-exon genes have `tx_pos` or `tx_neg` (not intergenic) and `gene_strand != 0` (single-strand). They have `n_spliced = 0`. The strand decomposition estimates their RNA fraction. If expressed at high levels, the decomposition shows high $\hat{f}_{\text{RNA}}$ and seeds them as RNA ($\gamma = 0$). They never enter the gDNA seed via the intergenic path (because they have `tx_pos` or `tx_neg`). They never enter via the density-eCDF path (which is deleted).

### Antisense overlap regions

Regions with both `tx_pos` and `tx_neg` have `gene_strand == 0`. Strand decomposition is not applied. These regions get $\gamma = 0.5$ (ambiguous). The EM classifies them using density and FL only. This is correct — with RNA on both strands, strand-based gDNA estimation is impossible.

---

## Files Changed

| File | Changes |
|------|---------|
| `src/rigel/native/resolve_context.h` | Add `counts_unique[]` array. Accumulate unique-mapper counts separately. |
| `src/rigel/native/bam_scanner.cpp` | Expose `counts_unique` in the result dictionary alongside `counts`. |
| `src/rigel/pipeline.py` | Extract `region_counts_unique` from scan result. Pass to `calibrate_gdna()`. Remove `density_percentile`, `min_gdna_regions` param forwarding. |
| `src/rigel/calibration.py` | Rewrite `_seed_initial_partition()` with strand decomposition + intergenic seeding. Modify `_m_step()` for length-weighted MAP updates. Delete algebraic fallback. Add `_strand_decompose()` helper. Remove `GDNA_INIT_*` constants. |
| `src/rigel/config.py` | `CalibrationConfig`: remove `density_percentile`, `min_gdna_regions`. Add `prior_weight_fraction: float = 0.01`. |
| `docs/calibration/CALIBRATION_ALGORITHM.md` | Update to reflect V2 algorithm. |

---

## Implementation Order

```
Phase 1: Length-weighted M-step          ← independent, do first. Low risk.
Phase 2: Strict region eligibility       ← C++ change (unique-only counts).
          (can parallel with Phase 1)       Requires recompilation.
Phase 3: Strand-decomposed seeding       ← depends on Phase 2 (uses unique counts).
Phase 4: Empirical Bayes MAP priors      ← depends on Phases 1 + 3 (uses seed + length-weights).
Phase 5: Delete fallback + thresholds    ← depends on Phase 4 (priors stabilize EM).
```

Phases 1 and 2 are independent and can be developed in parallel. Phases 3–5 are sequential.

---

## Verification Plan

1. **Unit tests** (`pytest tests/`): All ~1029 tests must pass. Golden outputs will change; regenerate with `--update-golden`.

2. **Scenario-aligned tests** (`pytest tests/scenarios_aligned/`): All 49 tests must pass. These exercise small-genome scenarios that previously triggered the fallback.

3. **Edge case validation** (manual debug scripts):
   - (a) Pure RNA / zero gDNA → $\gamma \approx 0$, $\lambda_G \approx 0$
   - (b) Pure gDNA / zero RNA → $\gamma \approx 1$
   - (c) 5-region mini genome → EM converges with Bayesian priors
   - (d) Unstranded library ($SS = 0.5$) → strand signal zeros out, no crashes
   - (e) Single-exon gene scenario → gene not misclassified as gDNA

4. **Synthetic benchmark sweep** (`scripts/benchmark/configs/*.yaml`): Compare gDNA relative error against V1 baseline. Focus on low-contamination scenarios.

5. **Full-scale Armis2 benchmark** (if available): Run all 13 simulation conditions. Compare calibration $\pi$, $\gamma$ distributions, and downstream mRNA quantification accuracy.

---

## Single New Hyperparameter

`prior_weight_fraction: float = 0.01` in `CalibrationConfig`.

This controls how strongly the empirical Bayes prior resists the data. At 0.01 (1%), the prior contributes the equivalent of 1% of the total eligible basepairs as virtual observations. Larger values make the prior more conservative (resist data more); smaller values make it weaker (closer to MLE).

All other thresholds (`GDNA_INIT_MIN_REGIONS`, `GDNA_INIT_DENSITY_PERCENTILE`, `SS > 0.55`) are deleted.

---

## Deferred Items

1. **Spatial smoothing.** The boundary-sweep algorithm cuts the genome at biological discontinuities (intron-exon boundaries). Smoothing across these boundaries would smear the sharp density transitions the EM exploits to separate RNA from gDNA. Spatial smoothing requires careful treatment of boundary types and hard constraints. Deferred until proven necessary by benchmarking.

2. **Parametric FL models.** Length-weighting may stabilize the non-parametric histogram FL model by reducing the influence of noisy small regions. If FL histograms remain unstable at low read depths after length-weighting, parametric models (log-normal) can be added as a separate change.

3. **Strand decomposition as an E-step signal.** The strand decomposition is used for *seeding* only. It could additionally be incorporated as a 4th LLR signal in the E-step (alongside density, strand, FL). However, the strand LLR already captures strand balance information via the Beta-Binomial model. Adding the decomposition as a separate E-step signal risks double-counting strand information. Deferred until ablation studies can evaluate the marginal benefit.

4. **Per-chromosome calibration.** The V2 algorithm is global (all chromosomes pooled). Sex chromosomes, mitochondria, and control sequences may have different contamination levels. A hierarchical model with per-chromosome gDNA densities is future work.
