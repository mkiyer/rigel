# Calibration Algorithm V3: Continuous gDNA–RNA Deconvolution

## 1. Core Principle

**Calibration's sole job is to separate genomic DNA (gDNA) from Total RNA.**

Calibration does not distinguish between mature RNA (mRNA) and nascent RNA (nRNA). That is the downstream EM's job — the per-locus EM solver has access to exact exon/intron geometries and per-fragment splice/overhang/mismatch likelihoods that calibration does not. Calibration operates at the coarser level of genomic regions and aggregate statistics.

The output of calibration is simple: for every genomic region, a continuous estimate of how many unspliced fragments are gDNA:

$$\mathbb{E}[N_{\text{gDNA},r}]$$

Not a probability. Not a binary classification. A physical count — the expected number of gDNA fragments in region $r$. This estimate flows into the per-locus EM as a Dirichlet prior, scaled by a discount factor.

---

## 2. Why gDNA vs Total RNA is the Right Boundary

The three fragment types in a BAM file are:

| Type | Origin | Strand | Splice | Density |
|------|--------|--------|--------|---------|
| **mRNA** | Processed transcript | Follows SS | Spliced at annotated junctions | Proportional to transcript expression |
| **nRNA** | Nascent transcript | Follows SS | Unspliced (intron-retained) | Proportional to transcription rate |
| **gDNA** | Genomic DNA contamination | Symmetric (50/50) | Unspliced | Uniform background rate $\lambda_G$ |

mRNA and nRNA share a critical property: **both follow the library's strand specificity (SS)**. gDNA does not — it is double-stranded and contributes equally to both strands. This makes gDNA vs Total RNA a clean separation:

- **gDNA**: strand-symmetric, density $\approx \lambda_G$ (uniform background)
- **Total RNA (mRNA + nRNA)**: strand-biased toward SS, density correlated with local gene expression

The downstream EM handles mRNA vs nRNA separation using exon/intron geometry, splice junction evidence, and fragment-length models — signals that require per-fragment resolution and transcript-level modeling, not the aggregate region-level view that calibration has.

---

## 3. Available Signals

### 3.1 Strand Balance (Primary Signal)

For stranded libraries ($SS > 0.5$), the strand ratio of unspliced reads in a gene-strand-annotated region directly reveals the gDNA fraction. This is an analytical computation, not a fitted parameter.

**Properties:**
- Degrades continuously as $SS \to 0.5$ (the $1/(SS - 0.5)$ amplification)
- Vanishes entirely at $SS = 0.5$ (unstranded)
- Requires `gene_strand != 0` (unambiguous strand annotation)
- No iteration needed — direct formula

### 3.2 Read Density (Secondary Signal)

gDNA contamination produces a roughly uniform background density $\lambda_G$ across the genome. Expressed genes have higher density. Unexpressed regions have density $\approx \lambda_G$.

**Properties:**
- Works for both stranded and unstranded libraries
- Requires identifying the background density level
- Noisier than strand (density is confounded by nRNA and variable gene expression)
- Provides the only gDNA signal for unstranded libraries

### 3.3 Fragment Length (Derived)

gDNA and RNA typically have different fragment-length (FL) distributions — gDNA reflects the library-prep sonication/tagmentation profile, while RNA FL is influenced by transcript structure. The gDNA FL model is not used for calibration itself; it is passed to the downstream scoring module as a per-fragment likelihood signal.

---

## 4. The Strand Pathway

### 4.1 Analytical Decomposition

For a region $r$ with unambiguous gene strand (`gene_strand` $\neq 0$), the unspliced reads are a mixture of gDNA and Total RNA:

- RNA: contributes fraction $SS$ to sense, $(1 - SS)$ to antisense
- gDNA: contributes $0.5$ to each strand

Let $N_r$ = total unspliced fragments, $N_{\text{anti},r}$ = antisense unspliced fragments (gene-relative). Then:

$$N_{\text{anti},r} = N_{\text{RNA},r} \cdot (1 - SS) + N_{\text{gDNA},r} \cdot 0.5$$

Solving:

$$\boxed{\mathbb{E}[N_{\text{gDNA},r}^{\text{strand}}] = \max\!\left(0,\; \frac{N_{\text{anti},r} - N_r \cdot (1 - SS)}{SS - 0.5}\right)}$$

The $\max(0, \cdot)$ clamp handles sampling noise that can produce negative estimates.

**Conversion from R1-position to gene-relative antisense:**

Under R1-antisense convention (dUTP/TruSeq Stranded):
- `gene_strand == +1`: $N_{\text{anti}} = N_{\text{pos}}$ (R1 maps to + strand = antisense for + gene)
- `gene_strand == -1`: $N_{\text{anti}} = N_{\text{neg}}$

This conversion already exists in `compute_sense_fraction()` ([calibration.py:280–312](../../src/rigel/calibration.py#L280)).

### 4.2 Noise Characteristics

Under Binomial sampling of strand assignment:

$$\text{Var}[\mathbb{E}[N_{\text{gDNA},r}^{\text{strand}}]] = \frac{N_r \cdot SS \cdot (1 - SS)}{(SS - 0.5)^2}$$

| SS | Noise amplification | Behavior |
|----|---------------------|----------|
| 1.00 | $\times 0$ | Perfect: antisense = gDNA/2 exactly |
| 0.95 | $\times 1.2$ | Excellent |
| 0.90 | $\times 2.3$ | Good |
| 0.75 | $\times 12.0$ | Moderate |
| 0.55 | $\times 990$ | Useless |
| 0.50 | $\div 0$ | Undefined |

For individual low-count regions, the per-region estimate is noisy. This is acceptable because: (a) the discount factor $\beta$ (Section 7) reduces the prior's influence, and (b) the per-locus EM has its own per-fragment evidence to correct errors.

### 4.3 Deriving Global $\lambda_G$ from Strand Estimates

When the strand pathway is available ($SS > 0.5$), we can derive the global gDNA background density directly:

$$\lambda_G^{\text{strand}} = \frac{\sum_r \mathbb{E}[N_{\text{gDNA},r}^{\text{strand}}]}{\sum_r L_r}$$

where the sum is over all eligible regions with `gene_strand` $\neq 0$. This is a strand-derived density estimate — no percentile, no mixture model. It aggregates the analytical decomposition across many regions, where the per-region noise averages out.

**This $\lambda_G$ is then used by the density pathway** (Section 5) for regions where strand decomposition is unavailable.

### 4.4 Regions Excluded from Strand Decomposition

| Condition | Reason |
|-----------|--------|
| `gene_strand == 0` (antisense overlap) | Sense vs antisense is undefined when transcripts exist on both strands |
| `gene_strand == 0` (intergenic) | No gene strand annotation at all |
| $SS \leq 0.5$ | No strand signal (unstranded library) |

These regions use the density pathway exclusively (Section 5).

---

## 5. The Density Pathway

### 5.1 Purpose

The density pathway estimates $\mathbb{E}[N_{\text{gDNA},r}]$ for regions where the strand decomposition is unavailable or unreliable. This includes:

1. All regions when the library is unstranded ($SS = 0.5$)
2. Regions with `gene_strand == 0` (antisense overlap or intergenic)

The core idea: gDNA produces a roughly uniform background density $\lambda_G$ across the genome. If a region's unspliced density exceeds $\lambda_G$, the excess is RNA. If it is at or below $\lambda_G$, the fragments are predominantly gDNA.

### 5.2 Estimating $\lambda_G$

**Case A: Stranded library ($SS > 0.5$).** Use the strand-derived $\lambda_G^{\text{strand}}$ from Section 4.3. This avoids circular reasoning — the density pathway uses a $\lambda_G$ anchored by the strand signal, not by density itself.

**Case B: Unstranded library ($SS = 0.5$).** The strand pathway produces no signal. We estimate $\lambda_G$ from the density distribution directly.

gDNA density is approximately uniform across the genome. RNA expression is highly skewed — a few genes dominate, while vast swaths of the genome are transcriptionally silent. The bottom quantile of the density distribution therefore represents regions closest to pure gDNA background.

$$\lambda_G^{\text{density}} = \text{percentile}_K\!\left(\left\{\frac{N_{\text{unspliced},r}}{L_r} : r \in \text{eligible}\right\}\right)$$

where $K$ is a low percentile (default: 10th). This uses **unspliced density** specifically, because spliced fragments are definitively RNA and would inflate the density of expressed regions.

**Robustness notes:**
- Length-weight the percentile computation: large regions provide more reliable density estimates. A 50,000bp region with 5 reads ($d = 10^{-4}$) is a more informative anchor than a 50bp region with 0 reads ($d = 0$).
- If the library is very clean (near-zero gDNA), the bottom percentile may still reflect low-expression RNA rather than gDNA. This results in a mild overestimate of $\lambda_G$, producing slightly elevated gDNA priors. The per-locus EM will correct this — the per-fragment FL and geometry evidence override weak priors.

### 5.3 Per-Region gDNA Estimate from Density

$$\boxed{\mathbb{E}[N_{\text{gDNA},r}^{\text{density}}] = \min\!\left(N_{\text{unspliced},r},\; \lambda_G \cdot L_r\right)}$$

This says: the expected gDNA count in region $r$ is at most $\lambda_G \cdot L_r$ (the background rate times the region length), but cannot exceed the observed unspliced count.

---

## 6. Blending

For each region, we blend the strand and density estimates using a continuous weight derived from $SS$:

$$w_{\text{strand}} = (2 \cdot SS - 1)^2 = \frac{(SS - 0.5)^2}{0.25}$$

| SS | $w_{\text{strand}}$ | Interpretation |
|----|---------------------|----------------|
| 1.00 | 1.00 | Pure strand pathway |
| 0.90 | 0.64 | Mostly strand |
| 0.75 | 0.25 | Mostly density |
| 0.60 | 0.04 | Almost entirely density |
| 0.50 | 0.00 | Pure density pathway |

For regions with `gene_strand == 0` (no strand decomposition possible), $w_{\text{strand}} = 0$ regardless of library $SS$.

The blended estimate:

$$\boxed{\mathbb{E}[N_{\text{gDNA},r}] = w_r \cdot \mathbb{E}[N_{\text{gDNA},r}^{\text{strand}}] + (1 - w_r) \cdot \mathbb{E}[N_{\text{gDNA},r}^{\text{density}}]}$$

where $w_r = w_{\text{strand}}$ when `gene_strand` $\neq 0$, and $w_r = 0$ otherwise.

---

## 7. gDNA Fragment-Length Model

The gDNA FL model enters the downstream fragment scoring as a per-fragment likelihood signal, helping the EM distinguish gDNA from RNA candidates with different insert sizes.

### 7.1 Source Fragments

**Stranded libraries ($SS > 0.5$):**

Antisense unspliced fragments from expressed single-strand regions are overwhelmingly gDNA. At $SS = 0.95$, approximately $\frac{0.5}{0.5 + 0.05} \approx 91\%$ of antisense unspliced fragments are gDNA. We collect fragment lengths from antisense unspliced unique-mapper fragments in regions with `gene_strand` $\neq 0$.

To maximize purity, optionally filter to regions where $\mathbb{E}[N_{\text{gDNA},r}] / N_r > 0.5$ (gDNA-majority regions), though for highly stranded libraries the antisense filter alone provides sufficient purity.

**Unstranded libraries ($SS = 0.5$):**

Without strand discrimination, we identify the highest-confidence gDNA fragments via density: fragments from the bottom $K$th percentile of density regions (the same regions used to estimate $\lambda_G$). These are the most likely gDNA regions. The FL distribution from these regions is a reasonable approximation, contaminated by some low-expression RNA.

**Fallback:** If neither approach yields sufficient observations (ESS < 50), fall back to the intergenic FL model trained during the BAM scan (the current `frag_length_models.intergenic`). This is the existing fallback.

### 7.2 Construction

The FL model is a discrete histogram (same as the current `FragmentLengthModel`):

1. Collect (fragment_length, weight) pairs from source fragments
2. Build a weighted histogram: `counts[fl] += weight`
3. Normalize and compute log-probabilities
4. Apply tail smoothing beyond the observed maximum

For the stranded pathway, the weight is 1.0 per antisense unspliced fragment (they're ≥90% gDNA, no per-fragment weighting needed). For the unstranded pathway, weight = 1.0 for fragments in low-density regions.

### 7.3 Required Data Change

The current `fl_table` stores `(region_id, frag_len)` for unique unspliced single-region fragments but does **not** store the fragment's strand. To build the gDNA FL model from antisense fragments specifically, we need to add `frag_strand` to the fl_table.

**C++ change in `RegionAccumulator::accumulate()`:** When appending to `fl_region_ids` / `fl_frag_lens`, also append the `frag_strand` (STRAND_POS or STRAND_NEG) to a new `fl_frag_strands` vector.

---

## 8. Aggregation to Loci

### 8.1 Per-Locus Expected gDNA

For each locus, we sum the per-region estimates across overlapping regions:

$$\mathbb{E}[N_{\text{gDNA},\ell}] = \sum_{r \in \text{overlap}(\ell)} \mathbb{E}[N_{\text{gDNA},r}]$$

This is a simple sum of physical counts — no weighted averaging, no γ fractions. If a locus spans 10 regions and 3 of them have nonzero expected gDNA, the locus inherits their sum.

The overlap query uses the existing `region_cr` cgranges interval tree, identical to the current `compute_gdna_locus_gammas()` mechanism.

### 8.2 Per-Locus Total Fragments

We also need the total fragment count for the locus to set the RNA prior budget:

$$N_{\text{total},\ell} = \sum_{r \in \text{overlap}(\ell)} N_{\text{total},r}$$

---

## 9. The EM Handoff

### 9.1 The Discount Factor

Calibration and the per-locus EM see overlapping evidence. Both observe the strand ratios: calibration uses aggregate strand statistics per region; the EM scores each fragment's strand individually. Using the full expected count as the Dirichlet prior would double-count the strand signal.

We apply a discount factor $\beta$ (default: 0.05) to convert calibration's physical count estimate into a prior of appropriate strength:

$$\alpha_{\text{gDNA}} = \beta \cdot \mathbb{E}[N_{\text{gDNA},\ell}]$$

**Interpretation:** The calibration estimate carries the evidential weight of $\beta \times N$ fragments. At $\beta = 0.05$, calibration's prior is worth 5% of the per-fragment evidence. The EM can easily override the prior when per-fragment likelihoods disagree, while still receiving a meaningful signal about the expected gDNA level.

**Scaling behavior:**

| Locus size | $\mathbb{E}[N_{\text{gDNA}}]$ | $\alpha_{\text{gDNA}}$ | Behavior |
|------------|-------------------------------|------------------------|----------|
| 10,000 frags, 10% gDNA | 1,000 | 50.0 | Strong gDNA prior — component competes robustly |
| 1,000 frags, 10% gDNA | 100 | 5.0 | Moderate prior — EM data dominates |
| 100 frags, 1% gDNA | 1 | 0.05 | Tiny prior — VBEM naturally suppresses |
| 1,000 frags, 0% gDNA | 0 | 0.0 | Zero prior — component fully disabled |

### 9.2 RNA Prior Budget

The RNA side of the Dirichlet prior is set symmetrically:

$$\alpha_{\text{RNA,total}} = \beta \cdot \left(N_{\text{total},\ell} - \mathbb{E}[N_{\text{gDNA},\ell}]\right)$$

This is distributed across individual RNA components (mRNA + nRNA per transcript) proportionally to their coverage-weighted representation, using the existing logic in `compute_ovr_prior_and_warm_start()`.

### 9.3 Warm Start

The gDNA component's initial $\theta$ is set proportionally:

$$\theta_{\text{gDNA}}^{(0)} = \frac{\alpha_{\text{gDNA}}}{\alpha_{\text{gDNA}} + \alpha_{\text{RNA,total}}} \cdot \theta_{\text{total}}$$

This replaces the current `locus_gamma / (1 - locus_gamma) * others` formula with an equivalent ratio derived from the physical count priors.

### 9.4 Interface Change

The C++ function `compute_ovr_prior_and_warm_start()` currently receives `locus_gamma` (float ∈ [0,1]) and `total_pseudocount` (float, default 1.0). Replace these with:

- `alpha_gdna` (float ≥ 0): the discounted expected gDNA count
- `alpha_rna_total` (float ≥ 0): the discounted expected RNA count

The Python side computes these from $\beta$ and the per-locus aggregates. The C++ side receives ready-to-use Dirichlet pseudo-counts.

---

## 10. EM Solver Changes

### 10.1 Delete the Binary Gate

Current code ([em_solver.cpp:1879–1881](../../src/rigel/native/em_solver.cpp#L1879)):

```cpp
if (locus_gamma == 0.0) {
    sub.prior[sub.gdna_idx] = 0.0;  // Disable gDNA
}
```

**Delete this block entirely.** When $\alpha_{\text{gDNA}} = 0$ (from $\mathbb{E}[N_{\text{gDNA}}] = 0$), the component's prior is already zero. The VBEM's `digamma(0)` returns $-\infty$, producing zero responsibility and zero posterior. The math disables the component without an `if` statement.

When $\alpha_{\text{gDNA}} > 0$ but very small (e.g., 0.001), the VBEM naturally suppresses it:

$$\psi(0.001) \approx -1000.4 \implies \text{log-weight} \ll \text{RNA components} \implies \text{responsibility} \approx 0$$

The component dies naturally. No gate needed.

### 10.2 Delete the VBEM Clamp Floor

Current code ([em_solver.cpp:946, 962](../../src/rigel/native/em_solver.cpp#L946)):

```cpp
double floor_i = std::max(prior[i], VBEM_CLAMP_FLOOR);
if (state_extrap[i] < floor_i) state_extrap[i] = floor_i;
```

The `VBEM_CLAMP_FLOOR = 0.1` prevents components from entering the "digamma absorbing regime." But this is the regime we **want** small components to be in. The entire point of VBEM is that digamma self-reinforces suppression of unlikely components.

**Replace with:** Clamp at `max(prior[i], EM_LOG_EPSILON)`. This prevents SQUAREM extrapolation from creating negative $\alpha$ values (invalid for the Dirichlet) while allowing arbitrarily small positive values.

```cpp
double floor_i = std::max(prior[i], EM_LOG_EPSILON);
if (state_extrap[i] < floor_i) state_extrap[i] = floor_i;
```

At $\alpha = 10^{-300}$, `digamma` returns $\approx -10^{300}$. This is a large negative number but perfectly finite. After max-subtraction in the E-step kernel ([em_solver.cpp:429](../../src/rigel/native/em_solver.cpp#L429)), the component's exp(log-weight) underflows to 0.0. The component receives zero responsibility and stays dead. Exactly correct.

### 10.3 Numerical Safeguards

The existing E-step kernel uses max-subtraction normalization before `exp()`:

```cpp
double max_val = ll[i * k] + log_weights[cidx[0]];
// ... find max ...
row[j] = val - max_val;  // max is 0, everything else is ≤ 0
```

This means `exp(row[j])` for a suppressed component is `exp(-1000)` ≈ 0.0 (IEEE double underflow to +0.0). No `NaN`, no `Inf`, no crash. The `fast_exp` SIMD paths handle this correctly via their early-zero-skip cutoff (`EXP_CUTOFF` at approximately −708).

No additional safeguards are needed. The existing numerical infrastructure handles arbitrarily small $\alpha$ values.

### 10.4 MAP-EM Compatibility

For MAP-EM (non-variational), the M-step uses:

$$\theta_k^{\text{MAP}} \propto \max(N_k + \alpha_k - 1, \; 0)$$

When $\alpha_k = 0.001$: $\theta_k \propto \max(N_k - 0.999, 0)$. A component needs at least $\sim 1$ fragment to survive. Natural and correct.

When $\alpha_k = 50$: $\theta_k \propto N_k + 49$. The prior adds 49 virtual fragments. Significant but subordinate to a locus with thousands of real fragments.

When $\alpha_k = 0$: $\theta_k \propto \max(N_k - 1, 0)$. The $-1$ from the Dirichlet mode formula means the component needs at least 1 fragment just to be nonzero. For a component with zero prior and zero data: $\theta_k = 0$.

Both MAP-EM and VBEM handle the physical-count prior correctly.

---

## 11. Data Flow

```
                     ┌─────────────────────────┐
                     │  BAM Scan (C++)          │
                     │  RegionAccumulator       │
                     │ ┌─────────────────────┐  │
                     │ │ counts[4 cols]       │──── per-region: n_unspliced_pos,
                     │ │                     │     n_unspliced_neg, n_spliced_pos,
                     │ │                     │     n_spliced_neg
                     │ ├─────────────────────┤  │
                     │ │ fl_region_ids       │──── unique unspliced single-region
                     │ │ fl_frag_lens        │     fragment lengths
                     │ │ fl_frag_strands [NEW]│──── strand per FL observation
                     │ └─────────────────────┘  │
                     │  StrandModel training    │──── SS (strand specificity)
                     │  FL model training       │──── RNA FL distribution
                     └─────────────────────────┘
                                │
                ┌───────────────┴───────────────┐
                ▼                               ▼
    ┌──────────────────────┐       ┌──────────────────────┐
    │ compute_region_stats │       │ FL table + strand    │
    │ n_pos, n_neg,        │       │ (per-fragment data)  │
    │ n_unspliced, n_spliced│      └──────────┬───────────┘
    │ gene_strand, L_r     │                  │
    └──────────┬───────────┘                  │
               │                              │
    ┌──────────┴───────────┐                  │
    │                      │                  │
    ▼                      ▼                  ▼
┌─────────────┐  ┌──────────────────┐  ┌─────────────────┐
│ Strand      │  │ Density Pathway  │  │ gDNA FL Model   │
│ Pathway     │  │                  │  │                 │
│ E[N_gDNA_r  │  │ λ_G estimate     │  │ from antisense  │
│  ^strand]   │  │ E[N_gDNA_r       │  │ unspliced frags │
│ per region  │  │  ^density]       │  │ (stranded) or   │
│             │  │ per region       │  │ low-density     │
│ → λ_G      │  │                  │  │ regions         │
│   ^strand   │  │                  │  │ (unstranded)    │
└──────┬──────┘  └────────┬─────────┘  └────────┬────────┘
       │                  │                     │
       └────────┬─────────┘                     │
                ▼                               │
       ┌───────────────────┐                    │
       │ Blending          │                    │
       │ w = (2·SS−1)²     │                    │
       │ E[N_gDNA_r]       │                    │
       └────────┬──────────┘                    │
                │                               │
                ▼                               ▼
    ┌──────────────────────┐       ┌────────────────────────┐
    │ Locus Aggregation    │       │ FragmentScorer         │
    │ Σ E[N_gDNA_r]        │       │ gdna_fl_log_prob       │
    │ per overlapping      │       │ (scoring.cpp)          │
    │ region               │       └────────────────────────┘
    └──────────┬───────────┘
               │
               ▼
    ┌──────────────────────┐
    │ Discount + Handoff   │
    │ α_gDNA = β · E[N]   │
    │ α_RNA  = β · (T−E[N])│
    └──────────┬───────────┘
               │
               ▼
    ┌──────────────────────┐
    │ C++ EM Solver        │
    │ No binary gate       │
    │ No VBEM clamp floor  │
    │ VBEM or MAP-EM       │
    │ Components die       │
    │ naturally via math   │
    └──────────────────────┘
```

---

## 12. Calibration Output

Replace the current `GDNACalibration` dataclass with a simpler structure:

| Field | Type | Purpose |
|-------|------|---------|
| `region_e_gdna` | `np.ndarray` (float64) | $\mathbb{E}[N_{\text{gDNA},r}]$ per region |
| `region_n_total` | `np.ndarray` (float64) | $N_{\text{total},r}$ per region (for RNA budget) |
| `gdna_fl_model` | `FragmentLengthModel` | gDNA fragment-length distribution |
| `lambda_gdna` | `float` | Global gDNA background density $\lambda_G$ |
| `strand_specificity` | `float` | Library SS (echoed for downstream reference) |

**Deleted fields:** `region_posteriors` (γ), `mixing_proportion` (π), `kappa_strand` (κ), `gdna_density_global`, `expressed_density`, `n_iterations`, `converged`, all diagnostic arrays.

No convergence fields because there is no iteration.

---

## 13. Algorithm Pseudocode

```python
def calibrate_gdna_v3(
    region_counts,     # DataFrame: n_unspliced_pos, n_unspliced_neg, n_spliced_pos, n_spliced_neg
    fl_table,          # DataFrame: region_id, frag_len, frag_strand [NEW]
    region_df,         # DataFrame: tx_pos, tx_neg, exon_pos, exon_neg, length
    strand_specificity,# float ∈ [0.5, 1.0]
    density_percentile=10,  # only used when SS = 0.5
) -> CalibrationResult:

    stats = compute_region_stats(region_counts, region_df)

    # Eligible: nonzero counts, nonzero length
    eligible = (stats["n_total"] > 0) & (stats["region_length"] > 0)

    # ── Strand Pathway ──
    has_strand = stats["gene_strand"] != 0
    N_unspliced = stats["n_unspliced"]
    N_anti = np.where(
        stats["gene_strand"] == 1, stats["n_pos"],
        np.where(stats["gene_strand"] == -1, stats["n_neg"], 0.0)
    )
    SS = strand_specificity

    e_gdna_strand = np.maximum(
        0.0,
        (N_anti - N_unspliced * (1 - SS)) / max(SS - 0.5, 1e-30)
    )
    e_gdna_strand[~has_strand | ~eligible] = 0.0

    # Global λ_G from strand decomposition
    if SS > 0.5:
        strand_mask = has_strand & eligible
        total_e_gdna = e_gdna_strand[strand_mask].sum()
        total_length = stats["region_length"][strand_mask].sum()
        lambda_G = total_e_gdna / max(total_length, 1.0)
    else:
        # Unstranded: percentile of unspliced density
        d_unspliced = N_unspliced[eligible] / stats["region_length"][eligible]
        lambda_G = np.percentile(d_unspliced, density_percentile)

    # ── Density Pathway ──
    e_gdna_density = np.minimum(
        N_unspliced,
        lambda_G * stats["region_length"]
    )
    e_gdna_density[~eligible] = 0.0

    # ── Blending ──
    w_strand = (2 * SS - 1) ** 2  # global weight
    w_per_region = np.where(has_strand, w_strand, 0.0)

    e_gdna = w_per_region * e_gdna_strand + (1 - w_per_region) * e_gdna_density

    # ── gDNA FL Model ──
    gdna_fl_model = build_gdna_fl_from_antisense(
        fl_table, stats, eligible, has_strand, SS
    )

    return CalibrationResult(
        region_e_gdna=e_gdna,
        region_n_total=stats["n_total"].astype(np.float64),
        gdna_fl_model=gdna_fl_model,
        lambda_gdna=lambda_G,
        strand_specificity=SS,
    )
```

```python
def compute_locus_priors(
    loci, index, calibration, beta=0.05
) -> tuple[np.ndarray, np.ndarray]:
    """Compute per-locus Dirichlet priors from calibration."""

    alpha_gdna = np.empty(len(loci), dtype=np.float64)
    alpha_rna  = np.empty(len(loci), dtype=np.float64)

    for li, locus in enumerate(loci):
        e_gdna_sum = 0.0
        n_total_sum = 0.0
        for ref, start, end in locus.merged_intervals:
            for _s, _e, rid in index.region_cr.overlap(ref, start, end):
                e_gdna_sum += calibration.region_e_gdna[rid]
                n_total_sum += calibration.region_n_total[rid]

        alpha_gdna[li] = beta * e_gdna_sum
        alpha_rna[li]  = beta * max(n_total_sum - e_gdna_sum, 0.0)

    return alpha_gdna, alpha_rna
```

---

## 14. Edge Cases

### Pure RNA (Zero gDNA)

Strand decomposition: $N_{\text{anti}} \approx N \cdot (1 - SS)$ everywhere. The formula yields $\mathbb{E}[N_{\text{gDNA}}] \approx 0$ per region. $\lambda_G \approx 0$. Density pathway: $\lambda_G \cdot L_r \approx 0$. Blended result: $\mathbb{E}[N_{\text{gDNA}}] \approx 0$ everywhere. EM prior $\alpha_{\text{gDNA}} \approx 0$. gDNA component naturally dies. Correct.

### Pure gDNA (Zero RNA)

No spliced reads. Unspliced reads are strand-symmetric: $N_{\text{anti}} \approx N/2$. Strand decomposition: $\mathbb{E}[N_{\text{gDNA}}] \approx N$ everywhere (since $N/2 - N \cdot (1-SS)/(SS-0.5) = N$). $\lambda_G$ is the global density. EM prior $\alpha_{\text{gDNA}} \gg \alpha_{\text{RNA}}$. gDNA component dominates. Correct.

### Unstranded Library ($SS = 0.5$)

$w_{\text{strand}} = 0$. Pure density pathway. $\lambda_G$ from percentile of unspliced density. $\mathbb{E}[N_{\text{gDNA},r}] = \min(N_r, \lambda_G \cdot L_r)$. Noisier than stranded but functional. The FL signal (Section 7) becomes the primary discriminator in the downstream EM.

### Single-Exon Genes

Have `tx_pos` or `tx_neg = True` (annotated), `n_spliced = 0`. With `gene_strand` $\neq 0$: strand decomposition applies. Expressed single-exon genes have strand-biased reads → $\mathbb{E}[N_{\text{gDNA}}]$ is small (correctly attributed to RNA). Unexpressed single-exon genes have strand-symmetric or zero reads → correctly attributed to gDNA.

### Antisense Overlap Regions

`gene_strand == 0`. Strand decomposition unavailable. These regions use density pathway only ($w_r = 0$). $\mathbb{E}[N_{\text{gDNA},r}] = \min(N_r, \lambda_G \cdot L_r)$.

### Exome Capture / Targeted Panels

No intergenic reads. **Not a problem.** The strand pathway operates on annotated regions (exonic and intronic within captured genes). The global $\lambda_G$ is derived from strand decomposition within captured regions, not from intergenic space. The density pathway applies $\lambda_G$ to regions without strand info. The entire algorithm works without a single intergenic read.

### Very Low gDNA Contamination

$\lambda_G \approx 0$. $\mathbb{E}[N_{\text{gDNA}}] \approx 0$ everywhere. $\alpha_{\text{gDNA}} \approx 0$. gDNA component suppressed in all loci. Correct — the tool doesn't hallucinate gDNA that isn't there.

### Very High gDNA Contamination

$\lambda_G$ is large. Most unspliced fragments attributed to gDNA. $\alpha_{\text{gDNA}}$ is large. gDNA component competes strongly in EM. The per-fragment FL and strand likelihoods still provide resolution — fragments with RNA-like FL and strand bias are correctly assigned to RNA despite the strong gDNA prior.

### Copy Number Aberrations

gDNA density is not perfectly uniform — amplified regions have higher gDNA density. The strand pathway handles this naturally (strand ratio is independent of copy number). The density pathway may underestimate gDNA in amplified regions (their density exceeds $\lambda_G$). This is a known limitation: the density pathway assumes uniform gDNA background.

---

## 15. What This Eliminates

| Removed | Reason |
|---------|--------|
| The calibration EM (50+ iterations, E-step, M-step, convergence check) | Replaced by analytical computation |
| Density Gaussians ($\mu_G, \sigma^2_G, \mu_E, \sigma^2_E$) | Not needed — gDNA density is a single scalar $\lambda_G$ |
| Beta-Binomial κ estimation | Not needed — strand signal is used analytically, not as a fitted distribution parameter |
| The algebraic fallback (70+ lines) | No edge cases requiring fallback |
| `GDNA_INIT_MIN_REGIONS = 100` | No minimum region count needed |
| `GDNA_INIT_DENSITY_PERCENTILE = 0.10` | Replaced by $\lambda_G$ from strand decomposition (stranded) or simple percentile (unstranded only) |
| `SS > 0.55` hard threshold | $w_{\text{strand}}$ degrades continuously |
| `VBEM_CLAMP_FLOOR = 0.1` | VBEM naturally suppresses small components |
| `if (locus_gamma == 0.0)` binary gate | Zero prior naturally disables component |
| `region_posteriors` (γ array) | Replaced by `region_e_gdna` (physical counts) |
| `mixing_proportion` (π) | Not needed — per-locus aggregation sums counts directly |
| `kappa_strand` | Never used downstream, always was diagnostic-only |
| `gdna_density_global`, `expressed_density` | Never used downstream |

---

## 16. New Parameters

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `beta` | 0.05 | Discount factor: prior strength = β × expected count. 0.05 means calibration evidence is worth 5% of per-fragment evidence. |
| `density_percentile` | 10 | Only used for unstranded libraries ($SS = 0.5$). Kth percentile of unspliced density for $\lambda_G$. |

Both parameters degrade gracefully:
- $\beta \to 0$: calibration has no influence, EM runs with flat priors (equivalent to no calibration)
- $\beta \to 1$: calibration evidence is fully counted (risk of double-counting, but no crash)
- `density_percentile` $\to 0$: $\lambda_G$ = minimum density (conservative)
- `density_percentile` $\to 50$: $\lambda_G$ = median density (aggressive gDNA estimation)

---

## 17. Implementation Plan

### Phase 1: C++ Data Collection (Low Risk)

**`src/rigel/native/resolve_context.h`:** Add `fl_frag_strands` vector to `RegionAccumulator`. Populate with `frag_strand` alongside `fl_frag_lens`. Expose in the result dictionary via `bam_scanner.cpp`.

**`src/rigel/pipeline.py`:** Extract `fl_frag_strands` from scan result. Add `frag_strand` column to `fl_table` DataFrame.

### Phase 2: Calibration Rewrite (Core Change)

**`src/rigel/calibration.py`:**
- Define new `CalibrationResult` dataclass (Section 12)
- Implement `calibrate_gdna_v3()` (Section 13) — the analytical strand + density pipeline
- Implement `build_gdna_fl_from_antisense()` — FL model from strand-selected fragments
- Delete the current EM machinery: `_e_step()`, `_m_step()`, `_seed_initial_partition()`, `estimate_kappa_marginal()`, iteration loop, convergence checks, algebraic fallback

**`src/rigel/config.py`:** Replace `CalibrationConfig` fields. Remove `max_iterations`, `convergence_tol`, `density_percentile`, `min_gdna_regions`, `min_fl_ess`. Add `beta: float = 0.05`, `density_percentile: int = 10`.

### Phase 3: Locus Prior Computation

**`src/rigel/locus.py`:** Replace `compute_gdna_locus_gammas()` with `compute_locus_priors()` (Section 13). Returns `(alpha_gdna, alpha_rna)` arrays instead of `locus_gammas`.

**`src/rigel/pipeline.py`:** Update the quantification phase to pass `alpha_gdna`, `alpha_rna` instead of `locus_gammas`.

### Phase 4: EM Solver Changes

**`src/rigel/native/em_solver.cpp`:**
- Delete `VBEM_CLAMP_FLOOR = 0.1`. Replace both clamp sites with `max(prior[i], EM_LOG_EPSILON)`.
- Modify `compute_ovr_prior_and_warm_start()`: accept `alpha_gdna` and `alpha_rna_total` instead of `locus_gamma` and `total_pseudocount`.
- Delete the `if (locus_gamma == 0.0) sub.prior[gdna_idx] = 0.0` block.
- Delete gDNA warm-start override (replaced by prior-ratio warm-start).

**`src/rigel/estimator.py`:** Update the Python-side EM dispatch to pass the new prior format.

### Phase 5: Test Updates

- Update all calibration tests to match new API and output format
- Regenerate golden outputs (`pytest tests/ --update-golden`)
- Verify all 1029+ tests pass
- Add specific test cases:
  - Pure RNA → $\alpha_{\text{gDNA}} \approx 0$
  - Pure gDNA → $\alpha_{\text{gDNA}}$ dominates
  - Unstranded → density pathway produces reasonable estimates
  - VBEM with tiny α → component naturally suppressed (no clamp)
  - VBEM with zero α → component disabled (no gate)

### Implementation Order

```
Phase 1 (C++ data) ──────────── independent, do first
Phase 2 (calibration rewrite) ─ depends on Phase 1
Phase 3 (locus priors) ──────── depends on Phase 2
Phase 4 (EM solver) ─────────── independent of Phases 2-3, can parallel
Phase 5 (tests) ─────────────── depends on all above
```

Phases 1 and 4 are independent and can be developed in parallel.

---

## 18. Verification

### Correctness Checks

1. **Conservation:** $\sum_r \mathbb{E}[N_{\text{gDNA},r}] \leq \sum_r N_{\text{unspliced},r}$. We never attribute more fragments to gDNA than exist.

2. **Monotonicity:** Higher global gDNA contamination → higher $\lambda_G$ → higher per-region $\mathbb{E}[N_{\text{gDNA}}]$ → higher locus $\alpha_{\text{gDNA}}$.

3. **SS continuity:** As $SS$ varies from 0.5 to 1.0, the algorithm output changes continuously. No discontinuities, no threshold effects.

4. **Zero gDNA invariant:** In a pure RNA sample, $\alpha_{\text{gDNA}} \approx 0$ for all loci.

### Benchmark Validation

1. **Synthetic sweep** (`scripts/benchmark/configs/*.yaml`): Run all configurations. Compare gDNA relative error against V1 baseline. Focus on low-contamination scenarios.

2. **Full-scale simulation** (Armis2): Run all 13 conditions. Compare calibration $\lambda_G$ against true gDNA density. Compare per-locus $\alpha_{\text{gDNA}}$ against oracle gDNA counts.

3. **VBEM stability test:** Construct loci with $\alpha_{\text{gDNA}} \in \{0, 10^{-6}, 10^{-3}, 1, 100, 10000\}$. Verify SQUAREM converges without the clamp floor. Verify small-α components are suppressed. Verify large-α components persist.

---

## 19. Summary

The V3 calibration algorithm replaces a 1400-line iterative EM mixture model with a ~100-line analytical computation. The key insights:

1. **Strand balance directly reveals gDNA.** No mixture model needed — the formula $\mathbb{E}[N_{\text{gDNA}}] = (N_{\text{anti}} - N \cdot (1-SS)) / (SS - 0.5)$ is exact given the strand model.

2. **Physical counts are the natural prior.** Dirichlet $\alpha$ values ARE pseudo-counts. Expressing the prior as discounted expected fragment counts gives it the right units, the right scale, and the right behavior across locus sizes.

3. **The VBEM math already handles component suppression.** The clamp floor and binary gate were workarounds for a parameterization problem (γ fraction × C=1 pseudocount). With physical-count priors, the math works as Bayesian theory intends: large α → component lives, small α → component dies.

4. **No intergenic dependency.** The strand decomposition operates on annotated regions (exonic, intronic). gDNA signals come from the strand balance within expressed genes, not from unexpressed genomic deserts. This works for exome capture, targeted panels, and any protocol that depletes intergenic space.

5. **Continuous degradation at low SS.** The blending weight $(2 \cdot SS - 1)^2$ smoothly transitions from pure strand (SS=1) to pure density (SS=0.5). No hard thresholds, no special cases.
