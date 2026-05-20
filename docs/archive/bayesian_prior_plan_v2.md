# Bayesian EM Prior Redesign Plan (v2)

**Date:** 2026-05-08
**Status:** Draft — Architecture and Theory Synthesis

---

## 1. The Theory: Biological Asymmetry

The fundamental theoretical breakthrough in this redesign is recognizing the biological asymmetry between the mixture model's components: **RNA is a biological signal**, while **gDNA is a physical artefact**.

In previous iterations, the prior system forced symmetry: it attempted to predict a total gDNA mass ($N_{gdna}$), converted it into a fraction ($\pi_{gdna} = N_{gdna} / N_{obs}$), and then allocated generic Dirichlet strength ($\alpha$) proportionally to both gDNA and RNA components via constants like `c_base`.

This was mathematically and biologically flawed:
1. **RNA is a regulated signal:** Transcription is locus-specific and dynamic. Global or regional levels do not predict whether a specific transcript is turned "on" or "off." An informative prior on RNA abundance is biologically presumptive and inappropriate. It must be driven entirely by locus-specific read data (the Likelihood).
2. **gDNA is a physical background:** Contamination from incomplete DNase digestion or capture probe off-target effects represents a physical sampling of the genome. It is subject to global averages and regional copy-number variations. Therefore, locoregional and global densities *are* highly predictive of exonic gDNA contamination.
3. **Double-counting evidence:** The previous locoregional prior used boundary reads inside the target locus to build the prior, then passed those exact same boundary reads into the EM, mechanically over-assigning reads to gDNA.

### The New Paradigm
We shift to an explicitly asymmetric, posterior predictive prior model:
*   We use independent background evidence to estimate the local gDNA density $\rho_l$.
*   We project that density onto the ambiguous exonic space to calculate an expected count of gDNA fragments: $\eta_g$.
*   We pass $\eta_g$ directly as the pseudo-count prior for gDNA.
*   We set an uninformative prior (0) for RNA, letting the Likelihood function and the in-locus data determine RNA expression.

```text
alpha_gdna = eta_g   (Informed by the physical background artefact model)
alpha_rna  = 0       (Uninformed, let the biological likelihood speak)
```

## 2. Design Principle: Strict Data Partitioning

To avoid double-counting, we must enforce conditional independence by formally splitting the observations into two sets:
*   **$C_l$ (Calibration Data):** External evidence used to estimate background density (e.g., global density, or flanking intergenic regions).
*   **$A_l$ (Ambiguous Target Data):** The reads routed to the EM locus to be deconvolved (including in-locus boundary and exonic fragments).

The prior must be derived *exclusively* from $C_l$.
The EM likelihood evaluates *only* $A_l$.

## 3. Solver Semantics: True Pseudocounts

Rigel's native solver handles Python `alpha` values as **physical pseudocount increments**.
*   **MAP-EM:** `alpha` acts as an extra weight equivalent to $\alpha - 1$.
*   **VBEM:** `alpha` acts as an increment on top of the Jeffreys baseline.

Therefore, we can directly pass the expected external count $\eta_g$ as `alpha_gdna` without translation into fractions or heuristics.

### Native Guardrail
Currently, the native engine uses `alpha_gdna <= 0` to disable the gDNA component. Since we want an uninformative RNA prior (`alpha_rna = 0`) and potentially a pure zero prior for gDNA when background is absent, we must decouple component eligibility from the prior mass scalar. 

## 4. Phased Implementation Plan

### Phase 0: Native Semantics Guardrail
- Decouple gDNA component eligibility from `alpha_gdna > 0` in the batch EM native path.
- Add a explicit boolean array input (e.g. `locus_enable_gdna: bool[n_loci]`) based on whether finite gDNA likelihood candidates exist.
- Ensure the gDNA warm-start behavior defaults to the coverage-derived guess when `alpha_rna == 0`, rather than zeroing out.

### Phase 1: Diagnostic and Schema Cleanup
- Add `gdna_prior_count` and `rna_prior_count` directly to the `PriorTable` and schema.
- Drop `c_base` and related heuristic tuning parameters from the schema and CLI arguments (mark them as deprecated/no-ops).
- Move `pi_gdna` inside Python to be an exclusively diagnostic, post-hoc metric, removing it entirely from the prior generation logic.

### Phase 2: Global-Only Independent Prior (The Baseline)
- Build the initial prior formulation using **global calibration densities** only. This is the simplest mathematically clean foundation.
- For each EM core locus:
  ```text
  eta_g = rho_global_intergenic * L_core_intergenic
        + rho_global_intron     * L_core_intron
        + rho_global_boundary   * B_core_boundary
        + rho_global_boundary   * L_core_exon_contained
  ```
- Set `alpha_gdna = eta_g` and `alpha_rna = 0`.
- Do not let in-locus boundary flux observations ($u_{left}, u_{right}$) contribute to this `eta_g`.

### Phase 3: Tests and Golden Update
- Update tests that previously expected `alpha_gdna + alpha_rna == c_loco`.
- Verify the zero-prior native eligibility.
- Ensure golden datasets reflect the new uninformative `alpha_rna` and pure pseudocount interpretation.

### Phase 4: Independent Flank Source (Production Destination)
- Implement `gdna_prior_source = independent_flank`.
- Compute physical background density exclusively from a set of flanking windows (`D_flank`) explicitly excluding the EM locus boundaries.
- Project this local independent density into `eta_g`, providing a capture-aware, CNV-aware local prior that still perfectly avoids double-counting $A_l$.
