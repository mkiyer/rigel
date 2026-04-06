# gDNA Calibration Algorithm

## Overview

The gDNA calibration module (`src/rigel/calibration.py`) classifies genomic regions into two categories—**expressed** (contains RNA) and **not expressed** (gDNA-only)—using a two-component mixture model. This classification produces per-region posteriors $\gamma_r \in [0,1]$ representing $P(\text{not expressed} \mid \text{data})$, which flow downstream into the locus-level EM solver as gDNA priors.

### Inputs

| Input | Description |
|-------|-------------|
| `region_counts` | Per-region fragment counts: `n_unspliced_pos`, `n_unspliced_neg`, `n_spliced_pos`, `n_spliced_neg` |
| `fl_table` | Fragment-length observations with `region_id` and `frag_len` |
| `region_df` | Region metadata: `tx_pos`, `tx_neg`, `length`, `ref`, optional `exon_pos`, `exon_neg` |
| `strand_specificity` | Library strand bias $SS \in [0.5, 1.0]$ learned from unique mappers |

### Outputs (GDNACalibration dataclass)

| Output | Description |
|--------|-------------|
| `region_posteriors` | Per-region $\gamma_r = P(\text{not expressed} \mid \text{data})$ |
| `gdna_density_global` | Global gDNA density (fragments/bp) |
| `kappa_strand` | Beta-Binomial strand concentration $\kappa$ |
| `gdna_fl_model` | Fragment-length distribution for gDNA |
| `mixing_proportion` | Global $\pi = P(\text{not expressed})$ |

---

## Regions

Regions are atomic genomic intervals produced by `build_region_table()` in `index.py`. The boundary-sweep algorithm projects all transcript exon boundaries (exon starts, exon ends, transcript starts, transcript ends) onto the genome and slices it into non-overlapping bins. Each bin has a fixed annotation context (which transcripts overlap, which are exonic vs intronic, which strands).

Region metadata includes:
- `tx_pos` / `tx_neg`: whether any plus/minus-strand transcripts overlap
- `exon_pos` / `exon_neg`: whether the region overlaps exonic sequence on each strand
- `length`: region span in bp

---

## Three Convergent Signals

The calibration EM uses three independent signals to distinguish gDNA from RNA:

### 1. Density (log-scale Gaussian)

gDNA is low-abundance background with roughly uniform genomic coverage. Expressed regions have higher fragment density. The model fits two Gaussians on $\log(\text{density})$:

$$\text{LLR}_{\text{density}} = \log \mathcal{N}(x \mid \mu_G, \sigma^2_G) - \log \mathcal{N}(x \mid \mu_R, \sigma^2_R)$$

### 2. Strand balance (Beta-Binomial)

gDNA produces reads equally on both strands ($p \approx 0.5$). RNA follows the library strand specificity ($p \approx SS$). The model uses gene-strand-normalized sense fractions:

- For `gene_strand == +1` (R1-antisense convention): $\text{sense\_frac} = 1 - \text{strand\_ratio}$
- For `gene_strand == -1`: $\text{sense\_frac} = \text{strand\_ratio}$

$$\text{LLR}_{\text{strand}} = \log \text{BetaBin}(k \mid n, \kappa/2, \kappa/2) - \log \text{BetaBin}(k \mid n, \kappa \cdot SS, \kappa \cdot (1-SS))$$

When $SS = 0.5$ (unstranded), the strand LLR is identically zero—no information.

### 3. Fragment-length shape

gDNA fragment lengths differ from RNA fragment lengths (gDNA tends to be longer with wider distribution). Uses shape-normalized histograms with Dirichlet smoothing to prevent total weight from biasing classification:

$$\text{LLR}_{\text{FL}} = \sum_b n_b \cdot \log(f_G(b) / f_R(b))$$

---

## Main EM Path

### Phase 1: Summary Statistics

`compute_region_stats()` produces per-region: `n_total`, `n_spliced`, `n_unspliced`, `density`, `strand_ratio`, `gene_strand`.

### Phase 2: Preprocessing

- Log-density: $\log(n/L + \epsilon)$ where $\epsilon = 1/\text{median}(L)$
- Sense fraction: gene-strand-normalized (see above)

### Phase 3: Initialization (`_seed_initial_partition`)

Two-phase seed partition:

**Expressed seed (hard constraint):** Any region with `n_spliced > 0` gets $\gamma = 0$. Splicing is definitive RNA evidence—this constraint is never relaxed during EM.

**gDNA seed (soft):** Among unspliced-only regions:
1. Build empirical eCDF of log-density from spliced (expressed) regions
2. Compute quantile rank for each unspliced region against the expressed density distribution
3. Regions below `density_percentile` (default 10th percentile) become gDNA candidates
4. Guarantee minimum seed size: if fewer than `GDNA_INIT_MIN_REGIONS` (100), take the bottom-N unspliced by quantile

**Supplementary strand filter:** Regions with symmetric strand balance (sense_frac $\approx 0.5 \pm 0.1$ when gene_strand $\neq 0$) are additionally seeded as gDNA.

Initial Gaussian parameters ($\mu$, $\sigma^2$) are computed from the seed weights. Initial FL models are built weighted by $\gamma$. Shared $\kappa$ is estimated via `estimate_kappa_marginal`.

### Phase 4: EM Iteration (max 50 iterations)

```
for iteration in range(max_iterations):
    # E-step: compute posterior γ from density + strand + FL LLRs
    log_odds = log_prior_odds + LLR_density + LLR_strand + LLR_FL
    γ[soft] = sigmoid(log_odds)
    γ[spliced] = 0.0  # hard constraint preserved
    
    # M-step: update π, density Gaussians, κ, FL models
    π_soft = mean(γ[soft regions])  # conditioned on unspliced
    update (μ_G, σ²_G), (μ_R, σ²_R) from weighted log-densities
    κ = estimate_kappa_marginal(...)
    gdna_fl, rna_fl = rebuild FL models weighted by γ
    
    # Convergence: |Δπ_soft| < 1e-4
    if converged: break
```

Key design: **π_soft** is the mixing proportion conditioned on unspliced-only regions (not global). Hard-constraint spliced regions ($\gamma = 0$) are excluded from the prior computation. This prevents spliced regions from diluting the prior for ambiguous regions.

### Phase 5: Final Outputs

One final E-step produces `final_gamma`. The gDNA FL model falls back to the intergenic FL model if the effective sample size is too small.

---

## Algebraic Fallback

### When It Triggers

```python
n_eligible = int(eligible.sum())
if n_eligible < min_gdna_regions:  # default GDNA_INIT_MIN_REGIONS = 100
```

The fallback fires when there are fewer than 100 eligible regions—too few for the EM density/strand/FL models to converge reliably. This commonly occurs in:
- Small genomes (test scenarios with 12–20 kb genomes produce only 4–6 regions)
- Targeted sequencing panels with limited genomic coverage
- Very sparse annotations

### What the Fallback Does

The fallback uses simple hard-constraint rules instead of iterative EM:

1. **Initialize** all $\gamma$ to 0.5 (neutral prior for regions without data)
2. **Spliced regions** ($n_{\text{spliced}} > 0$): $\gamma = 0$ (RNA)
3. **Unspliced-only regions** ($n_{\text{spliced}} = 0$, $n_{\text{total}} > 0$): classification depends on library strandedness:
   - **Stranded library** ($SS > 0.55$): Default to $\gamma = 0$ (RNA). Only reclassify to $\gamma = 1$ (gDNA) if a per-region binomial strand LLR test favors gDNA (balanced strand pattern with $\geq 5$ reads).
   - **Unstranded library** ($SS \leq 0.55$): Default to $\gamma = 1$ (gDNA), since no strand signal can distinguish RNA from gDNA.
4. **Strand LLR test** (stranded libraries only): For unspliced-only regions with gene-strand context and $\geq 5$ reads:

$$\text{LLR} = k \cdot \log\frac{SS}{0.5} + (n-k) \cdot \log\frac{1-SS}{0.5}$$

where $k$ = sense-strand count, $n$ = total unspliced count. If $\text{LLR} > 0$ (strand bias consistent with RNA), keep $\gamma = 0$. If $\text{LLR} \leq 0$ (balanced strands, consistent with gDNA), set $\gamma = 1$.

5. **Compute** $\pi_{\text{fallback}} = \text{mean}(\gamma_r \mid n_r > 0)$
6. **Compute** gDNA density from $\gamma > 0$ regions: $\lambda_G = \sum(\gamma_r \cdot n_r) / \sum(\gamma_r \cdot L_r)$

### Why the Default Matters

When $\gamma_{\text{locus}} > 0$ (even very small values like 0.0005), the gDNA component is marked **eligible** in the EM solver and participates in fragment assignment. The EM can then amplify even tiny gDNA priors into significant fragment counts if the gDNA scoring model is competitive with mRNA for some fragments. Setting the default to $\gamma = 0$ for stranded libraries prevents this: when all regions have $\gamma = 0$, the downstream `compute_gdna_locus_gammas` produces $\gamma_{\text{locus}} = 0$, which disables the gDNA component entirely.

### Previous Bug (pre-v0.4.0)

Before v0.4.0, the fallback unconditionally classified **all** unspliced-only regions as gDNA ($\gamma = 1$), regardless of library strandedness or strand balance.

**Why this was wrong:** For stranded RNA-seq ($SS \gg 0.5$), many unspliced-only regions contain RNA signal (nascent RNA, unspliced pre-mRNA, single-exon transcripts). These regions have strong strand bias ($\text{sense\_frac} \approx SS$), which is inconsistent with gDNA ($\text{sense\_frac} \approx 0.5$). Blindly labeling them as gDNA caused:

1. **Inflated $\pi$**: With both spliced ($\gamma = 0$) and unspliced ($\gamma = 1$) regions, $\pi_{\text{fallback}} \approx 0.5$, suggesting 50% of the genome is unexpressed—a gross overestimate for RNA-rich samples.

2. **gDNA component activation**: Even a single region with $\gamma = 1$ and one fragment could produce a non-zero $\gamma_{\text{locus}}$, enabling the gDNA component in the EM. Once enabled, the EM's gDNA component can compete for fragments and grow far beyond the calibration prior (see "Why the Default Matters" above).

3. **mRNA under-quantification**: In test scenarios with small genomes, this caused up to 36% of fragments to be incorrectly assigned to gDNA, producing severe mRNA under-estimation.

**The fix:** For stranded libraries ($SS > 0.55$), the fallback now defaults unspliced-only regions to RNA ($\gamma = 0$) instead of gDNA. Only regions with enough reads ($\geq 5$) and balanced strand patterns (LLR $\leq 0$) are classified as gDNA. This ensures the gDNA component is disabled when there is no genuine gDNA evidence, while preserving correct classification for unstranded libraries.

---

## Edge Cases and Known Limitations

### Pure RNA (no gDNA contamination)

When there is zero gDNA, the calibration should produce $\gamma \approx 0$ everywhere and $\lambda_G \approx 0$. The main EM handles this naturally: with no low-density symmetric regions, the gDNA component gets negligible weight. The fallback handles this via the strand-aware refinement (stranded regions → $\gamma = 0$).

### Pure gDNA (no RNA)

When there is no RNA expression, the calibration should produce $\gamma \approx 1$ everywhere. Edge case: no spliced regions exist, so `_seed_initial_partition` triggers the degenerate path ("Fewer than 2 spliced regions found. Assuming pure gDNA..."). The M-step safety net (`w_e_sum < 1.0`) handles the empty RNA component.

### Small genomes / few regions

The `GDNA_INIT_MIN_REGIONS = 100` threshold is conservative. With <100 regions, the Gaussian density model is underpowered, so the fallback is preferred. Future work could lower this threshold or use a non-parametric density model.

### Unstranded libraries ($SS \approx 0.5$)

The strand signal vanishes. The fallback's strand refinement is disabled ($SS \leq 0.55$). Classification relies solely on splicing. Unspliced-only regions default to $\gamma = 1$, which is a reasonable assumption for unstranded data where strand cannot distinguish RNA from gDNA.

---

## Downstream Flow

### Per-locus gDNA posteriors (`compute_gdna_locus_gammas` in locus.py)

For each locus, an aggregate $\gamma_{\text{locus}}$ is computed as the fragment-weighted average of overlapping region posteriors:

$$\gamma_{\text{locus}} = \frac{\sum_{r \in \text{overlap}} \gamma_r \cdot n_r}{\sum_{r \in \text{overlap}} n_r}$$

Falls back to `calibration.mixing_proportion` if no regions overlap.

### EM prior construction (`compute_ovr_prior_and_warm_start` in em_solver.cpp)

The locus-level $\gamma$ sets the gDNA component's prior budget. A locus with $\gamma \approx 0$ gets near-zero gDNA prior; one with $\gamma \approx 1$ gets a strong gDNA prior. This is the mechanism by which calibration controls fragment assignment.

---

## Open Issues

1. **Algebraic fallback is still crude**: Even with the strand fix, unspliced-only regions in unstranded libraries default to gDNA. A more sophisticated fallback could use density quantiles or fragment-length shape.

2. **`GDNA_INIT_MIN_REGIONS` threshold**: The value of 100 was chosen conservatively. It may be possible to run EM reliably with fewer regions using regularization or informative priors.

3. **Fragment-length deconvolution**: The FL signal is the weakest discriminator and can be noisy with small sample sizes. The shape normalization helps but doesn't eliminate the issue.

4. **Region granularity**: The boundary-sweep algorithm can produce very small regions at gene boundaries. These have high count variance and can be unreliable for density estimation.

5. **Multi-chromosome consistency**: Calibration is global (all chromosomes pooled). Chromosome-specific contamination patterns are not modeled.

---

## Key Constants

| Constant | Value | Location | Purpose |
|----------|-------|----------|---------|
| `GDNA_INIT_MIN_REGIONS` | 100 | `calibration.py:67` | Minimum regions for EM; below triggers fallback |
| `GDNA_INIT_DENSITY_PERCENTILE` | 0.10 | `calibration.py:66` | Bottom density quantile for gDNA seed |
| `_EPS` | 1e-12 | `calibration.py:62` | Numerical stability constant |
| `max_iterations` | 50 | `calibrate_gdna()` default | EM iteration limit |
| `convergence_tol` | 1e-4 | `calibrate_gdna()` default | $\|\Delta\pi\|$ convergence criterion |
| `min_fl_ess` | 50 | `calibrate_gdna()` default | Minimum FL effective sample size |
| `kappa_strand` | 2.0 | fallback default | Fixed $\kappa$ when EM isn't run |
