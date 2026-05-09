# Strand-Aware SRD Calibration Plan

## 1. Problem Statement

Rigel's Bayesian EM relies on accurate quantification of the global gDNA density through Simple Regional Deconvolution (SRD). Recent removal of heuristic scaling limits (Phase 2 of the prior redesign) revealed a fundamental flaw: overlapping nascent RNA (nRNA) fragments in intronic regions are systematically misclassified as gDNA by the SRD scanner. 

This causes severe double-counting, inflating the global gDNA prior $C_g$ and destroying the EM assignments in high-nRNA synthetic scenarios. The goal is to mathematically isolate gDNA counts ($N_g$) from RNA counts ($N_r$) during the initial calibration scan using the library's **Strand Specificity (SS)**.

## 2. Theoretical Foundation

### 2.1 The Probabilistic Deconvolution

Let $N$ be the total observed fragments mapping to a completely unambiguous stranded region (i.e., exactly one host gene strand, with no antisense overlapping genes).
Let $SS \in [0.5, 1.0]$ be the library's measured Strand Specificity (the empiric probability that a captured RNA molecule correctly maps to the sense strand of its originating transcript).
Let $N_g$ be the unobserved total number of gDNA fragments in this region.
Let $N_r$ be the unobserved total number of RNA fragments in this region.

When gDNA is fragmented and sequenced, it produces fragments equally from both Watson and Crick strands. Therefore, the probability of a gDNA fragment mapping to the annotated "Sense" strand is exactly $0.5$.

We can write the expected counts for Sense ($N_S$) and Antisense ($N_A$) mapped fragments:

$$E[N_S] = N_r \cdot SS + N_g \cdot 0.5$$
$$E[N_A] = N_r \cdot (1 - SS) + N_g \cdot 0.5$$

We can solve this system to isolate $N_r$:

$$E[N_S] - E[N_A] = N_r \cdot (2 \cdot SS - 1)$$
$$N_r = \frac{N_S - N_A}{2 \cdot SS - 1}$$

Substituting $N_r$ back into the total mass ($N_S + N_A = N_r + N_g$) gives the exact gDNA mass:

$$N_g = (N_S + N_A) - \frac{N_S - N_A}{2 \cdot SS - 1}$$

### 2.2 Boundary Conditions & Robustness

*   **Perfect Stranding ($SS = 1.0$):**
    $$N_g = (N_S + N_A) - (N_S - N_A) = 2 \cdot N_A$$
    This matches biological intuition seamlessly: if RNA is perfectly stranded, then *all* antisense reads must originate from unstranded gDNA contamination. Since gDNA is 50/50 symmetric, the total regional gDNA mass is simply twice the antisense mass.

*   **Unstranded Data ($SS \to 0.5$):**
    As $SS$ approaches $0.5$, the denominator $(2 \cdot SS - 1)$ approaches zero, making the variables linearly dependent and the system unidentifiable. To make SRD mathematically robust across all datasets, we introduce a smooth **Strand Trust Weight ($W_{SS}$)** that gracefully degrades the estimator back to the naive $N_{total}$ baseline when the library is unstranded.

    $$W_{SS} = \text{clamp}\left(\frac{SS - 0.55}{0.35}, 0, 1\right)$$
    $$N_{g, \text{final}} = W_{SS} \cdot N_{g, \text{deconvolved}} + (1 - W_{SS}) \cdot N_{total}$$

## 3. Implementation Plan

The fix requires capturing the relative orientation between the mapped fragment and the genomic reference region, storing these counts in a separated payload buffer, and applying the deconvolution in the Python SRD arrays.

### Phase 1: Native C++ Payload Upgrade

1.  **Augment `RegionIndex`:** 
    *   Update `region_index.h/.cpp` to accept a `py::array_t<uint8_t> strands` array during initialization.
    *   Store a `std::vector<uint8_t> strands_` backing array mapping `RegionStrand` enum (`NONE=0`, `POS=1`, `NEG=2`, `AMBIG=3`).
    *   Add a `uint8_t strand(int32_t rid)` accessor.
2.  **Augment `CalibrationPayload`:** 
    *   In `calibration_types.h`, split the three core payload arrays into 3 orientation buckets (`sense`, `anti`, `ambig`).
    *   `per_region_counts_sense` / `_anti` / `_ambig`
    *   `u_left_sense` / `_anti` / `_ambig`
    *   `u_right_sense` / `_anti` / `_ambig`
3.  **Update `CalibrationAccumulator::observe`:**
    *   Extract the mapped strand of the fragment using the first exon block (`int32_t frag_strand = obs_exons[0].strand;`).
    *   Compare `frag_strand` to the region's `strand(rid)`.
    *   Route the count: `_sense` if they match, `_anti` if they are opposite, and `_ambig` if the region is unstranded (`NONE`) or bidirectionally transcribed (`AMBIG`).
4.  **Python Bindings:** 
    *   Expose the 9 new arrays through `native.py` and `_result.py` (`CalibrationScanPayload`).

### Phase 2: SRD Global Density Deconvolution (Python)

1.  **Deconvolution Matrix Math:** 
    *   Update `_arrays.py` to aggregate the sense/anti/ambig buckets.
    *   In `density_global.py`, update `_density_one_type(INTRON)` to implement the matrix math from Section 2.1.
2.  **Boundary Deconvolution:** 
    *   In `_density_exon_intron()`, apply the exact same deconvolution logic independently to `u_left` and `u_right` fluxes.
3.  **Fallback Preservation:** 
    *   The `_ambig` arrays (representing overlaps of sense/antisense transcripts, or unstranded intergenic space) mathematically cannot be deconvolved in bulk this way. They will always fall back to $N_g = N_{total}$.

### Phase 3: Validation and Regression Recovery

1.  **Test Execution:** Run `test_nrna_double_counting.py` in highly stranded modes (`s90`, `s100`). The total RNA error should drop from ~44% back to near zero.
2.  **Goldens Regeneration:** Regenerate golden threshold arrays for the test harness to codify the new expected Bayesian corrections.
