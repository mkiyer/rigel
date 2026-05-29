# Rigel Calibration Fusion & EM Integration Implementation Plan (v1)

This document establishes the step-by-step implementation playbooks and validation tests for integrating the **Density Model** with the **Strand-Deconvolution Model** (Phase 4), constructing the locus-level **Prior Allocation** table (Phase 5), and wiring the results into the C++ **MultiLocus EM solver** (Phase 6).

---

## PR 1: Probabilistic Fusion Engine (`src/rigel/calibration/integration.py`)

### 1.1 Objective & Rationale
We implement a zero-heuristic, unified posterior estimator that integrates:
- Regional opportunity/exposure parameters modeled as continuous Negative Binomial variables over intergenic or intronic backdrops ($P_{\text{dens}}(d)$).
- Local, alignment-pure strand distributions modeled as a BetaBinomial mixture of unstranded gDNA and stranded/unstranded transcript orientations ($L_{\text{strand}}(K_r \mid D_r = d, N_r)$).

By exactly evaluating the discrete log-posterior for region counts $N_r \le 200$, and leveraging a bounded Laplace approximation beyond that inflection point, we achieve perfect numerical stability, mathematical rigor, and high computational performance.

---

### 1.2 Mathematical Formulation

For each region $r$, let:
* $N_r = \text{total observed unspliced reads} = C_r + B_r$ (contained + boundary).
* $K_r = \text{observed sense-like folded count}$.
* $D_r = \text{latent gDNA count} \in [0, N_r]$.
* $R_r = N_r - D_r = \text{latent RNA count} \in [0, N_r]$.

The unified log-posterior surface for the gDNA count $D_r = d$ is:
$$
\ell(d) = \log P_{\text{dens}}(D_r = d) + \log L_{\text{strand}}(K_r \mid D_r = d, N_r)
$$

Where:
* **Density Prior:** $P_{\text{dens}}(d) = \text{NB}(d; r = \alpha_{\text{post}, r}, p = p_{\text{nb}, r})$ from the `DensityEvidence` structure containing Gamma-conjugate parameters.
* **Strand Likelihood:** $L_{\text{strand}}(K_r \mid d, N_r)$ represents the marginal probability of observing $K_r$ sense reads when $d$ counts are gDNA and $r = N_r - d$ counts are RNA:
$$
L_{\text{strand}}(K_r \mid d) = \sum_{j=\max(0, K_r - d)}^{\min(R_r, K_r)} \text{Binomial}(j; \text{trials}=R_r, p=p_{r1\_\text{sense}}) \cdot \text{BetaBinomial}(K_r - j; \text{trials}=d, \alpha=\frac{\kappa_d}{2}, \beta=\frac{\kappa_d}{2})
$$

---

### 1.3 Implementation Steps

#### Step 1: Define `FusedRegionGdnaEvidence` Struct
Create a dataclass in a new file `src/rigel/calibration/integration.py`:
```python
@dataclass(frozen=True)
class FusedRegionGdnaEvidence:
    mean_count: np.ndarray          # float32 array of shape (R,) - E[d_r | data]
    upper_count: np.ndarray         # float32 array of shape (R,) - 95% upper credible limit
    tail_probability: np.ndarray    # float32 array of shape (R,) - Probability(d_r >= expected_density)
    flags: np.ndarray               # uint8 array of shape (R,)
```

#### Step 2: Implement Discrete Posteriors ($N_r \le 200$)
Evaluate $\ell(d)$ for all $d \in \{0, 1, \dots, N_r\}$:
1. Compute the array of log-prior density probabilities: `nbinom.logpmf(d, alpha_post, p_nb)`.
2. Compute the log-strand-likelihood terms using the exact convolution defined in `_exact_posterior_R` under `strand_deconv.py`. Use the variable mapping $D_r = n - r$ to obtain the gDNA perspective.
3. Compute the normalized posterior probabilities: $P_k = \exp(\ell(k) - \ell_{\max}) / \sum_j \exp(\ell(j) - \ell_{\max})$.
4. Extract expected count $\mathbb{E}[D_r] = \sum_k k \cdot P_k$.
5. Find the upper $1 - \alpha$ confidence bound: the smallest $M \in [0, N_r]$ such that $\sum_{k=0}^{M} P_k \ge 0.95$.

#### Step 3: Implement Bounded Laplace Approximation ($N_r > 200$)
1. Find the modal value $d^* = \text{argmax}_{d \in [0, N_r]} \ell(d)$ using a 1D vectorized golden-section or bounded ternary search.
2. Estimate the negative Hessian $H$ (second derivative of $- \ell(d)$ evaluated at $d^*$) via centered finite differences:
   $$H \approx -\frac{\ell(d^* + \epsilon) - 2\ell(d^*) + \ell(d^* - \epsilon)}{\epsilon^2}$$
3. Define the truncated normal distribution approximation: $d \sim \mathcal{N}\left(d^*, \sigma^2 = \frac{1}{H}\right)$ bounded within $[0, N_r]$.
4. Retrieve the analytical truncated mean, standard deviation, and cumulative density function (CDF) to determine the $95\%$ upper quartile.

---

### 1.4 Verification & Testing Plan
Create `tests/test_integration.py` to assert correct deconvolution behavior across physical extremes:
- **Test Case 1 (Strand-Specific Pure RNA / High Contamination):** Under high strand specificity, a massive $100\%$ sense-strand region (e.g. $10,000$ sense, $0$ antisense with $p_{r1\_\text{sense}} = 0.99$) must shrink the estimated gDNA to $0$ regardless of areally predicted background density.
- **Test Case 2 (Unstranded Library Fallback):** For unstranded settings ($p_{r1\_\text{sense}} = 0.50$), the joint posterior must match $P_{\text{dens}}(d)$ exactly, reverting fully to the Gamma rate prior and intronic/intergenic background calibration context.

---
---

## PR 2: Footprint-Aware Prior Allocator (`src/rigel/calibration/prior.py`)

### 2.1 Objective & Rationale
Prevent **hotspot smearing** (where high density or transcript abundance on one strand is misattributed to overlapping or neighboring locus windows on the opposite strand) by calculating priors using fine-region native interval intersections rather than raw sequence length.

---

### 2.2 Mathematical Metrics

*   **gDNA Continuous Relative Exposure ($A_r$):**
    For each fine region $r$, we calculate continuous scaling local exposure factors relative to the library reference density ($\rho_{\text{ref}}$):
    $$A_r = \frac{\mu_{D_r}}{\rho_{\text{ref}} \cdot L_{\text{eff}, r}^{\text{raw}}}$$
    Where $L_{\text{eff}, r}^{\text{raw}}$ is the raw physical base-pair span of the region.

*   **Locus-level gDNA Prior Counts ($D_{\text{EM}}^l$):**
    We project local estimates onto overlapping EM loci footprints:
    $$D_{\text{EM}}^l = \sum_{r \in \text{footprint}(l)} \theta_{r \to l} \cdot \mu_{D_r}$$
    Where $\theta_{r \to l}$ is the overlap ratio computed dynamically, capturing how many coordinate fragments map directly inside locus $l$.

*   **Locus-level gDNA Effective Length ($L_{\text{EM}}^l$):**
    $$L_{gDNA}^l = \sum_{r \in \text{footprint}(l)} A_r \cdot L_{\text{eff}, r}^{\text{raw}}$$
    By multiplying the local base scale raw span $L_{\text{eff}, r}^{\text{raw}}$ by the continuous exposure factor $A_r$, we ensure custom capture targets or highly localized hyper-exposed motifs receive proportionally shorter or longer scaled effective genomic bounds.

---

### 2.3 Implementation Steps

#### Step 1: Create `PriorTable` Dataclass
```python
@dataclass(frozen=True)
class PriorTable:
    gdna_prior_count: np.ndarray    # float32 of shape (L_loci,) - D_EM per locus
    gdna_eff_len: np.ndarray        # float32 of shape (L_loci,) - L_gDNA per locus
    enable_gdna: np.ndarray         # bool-mask array of shape (L_loci,)
    exposure_factors: np.ndarray    # float32 scaling array of shape (R_regions,)
```

#### Step 2: Establish the Regional Prior Filters
Exclude untrustworthy regions from generating exposure contributions by checking:
*   `FLAG_INELIGIBLE` or `FLAG_NEAR_UNSTRANDED` on the regions.
*   Assign $A_r = 1.0$ (flat exposure) for any non-anchor or flagged regions to protect against extreme statistical divisions.

#### Step 3: Implement Intersection Mapping
Intersect the coordinates of non-overlapping physical Locus footprints (provided by the `MultiLocus` object) with the $1000\text{bp}$ bins of the Ledger using interval indexing:
1. Traverse sorted interval coordinates.
2. For each region, distribute $\mu_{D_r}$ proportionally across intersecting loci according to exact coordinate coverage.
3. Sum the continuous exposure-weighted lengths to compute $L_{gDNA}^l$.

---

### 2.4 Verification & Testing Plan
Create `tests/test_prior.py`:
- Mock a locus layout containing:
  - Locus A (unexpressed) and Locus B (housing a massive 50,000 count capture target hotspot).
- Assert that under footprint-aware boundary matching, the estimated prior of Locus A receives exactly $0$ counts, and Locus B absorbs $100\%$ of the localized prior weight.

---
---

## PR 3: End-to-End Pipeline Wiring (`src/rigel/pipeline.py`)

### 3.1 Objective & Rationale
Remove the final `NotImplementedError` inside `quant_from_buffer` to fully wire the Python calibration modules with the highly optimized C++ Equivalence Class EM Solver, allowing Rigel to perform end-to-end transcript quantification over real and simulated datasets.

---

### 3.2 Implementation Steps

#### Step 1: Modify `src/rigel/pipeline.py`
Locate the block around line 1060 and replace the short-circuit exception. Replace it with the complete structural conduit:

```python
    # 1. Fuse local density estimates with alignment strand deconvolution
    fused_evidence = fuse_density_and_strand(
        calibration.density_evidence,
        calibration.region_gdna,
        calibration.calibration_context
    )

    # 2. Build multi-loci structures
    multi_loci = build_multi_loci(em_data, index)
    _assign_locus_ids(estimator, multi_loci)

    # 3. Build Prior Tables
    prior_table = assemble_priors(
        multi_loci, fused_evidence, em_data, index
    )

    # 4. Partition and execute parallel EM
    partitions = pr.build_partitions(em_data, multi_loci, estimator)
    _run_locus_em_partitioned(
        estimator, partitions, multi_loci, index,
        gdna_prior_count_em=prior_table.gdna_prior_count,
        gdna_eff_len=prior_table.gdna_eff_len,
        enable_gdna=prior_table.enable_gdna,
        em_config=em_config or EMConfig(),
        annotations=annotations,
        emit_locus_stats=emit_locus_stats
    )
```

#### Step 2: Register Outputs in Result Object
Ensure `PriorTable` and output parameters are registered as fields inside the final `CalibrationResult` named tuple so they are correctly flushed to output parquet or feather files.

---

### 3.3 Verification & Testing Plan
- Run the full synthetic integration scenarios:
  ```bash
  conda activate rigel
  pytest tests/test_pipeline_wiring.py -v
  pytest tests/test_golden_output.py -v
  ```
- If model updates require it, update the ground truth regression footprints:
  ```bash
  pytest tests/ --update-golden
  ```
- Trigger the complete benchmarking sweep configurations using the simulated standard:
  ```bash
  python -m scripts.benchmarking status -c scripts/benchmarking/configs/default.yaml
  ```
