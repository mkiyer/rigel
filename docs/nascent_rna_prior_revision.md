In the new architecture, nascent RNA is a physically distinct molecule. The nascent fraction ($\eta$) is a property of the **nRNA**, not the transcript. Therefore, your priors must be **per-nRNA** (size $N$), not per-transcript (size $T$).

Because an nRNA span exactly defines the transcription start and end sites, **the nRNA *is* the TSS Group.** The old hierarchy (`Transcript -> TSS Group -> Locus -> Global`) gracefully collapses into a beautifully simple 3-tier hierarchy:

1. **nRNA Level (Base):** The local evidence for this specific pre-mRNA span.
2. **Locus-Strand Level:** The pooled evidence for all nRNAs on the same strand in this locus.
3. **Global Level:** The global baseline for the entire dataset.

Here is exactly how to structure this mathematically and how to rewrite the Python code in `estimator.py`.

### 1. The Evidence Gathering

To build the prior for an nRNA, we need to compare its "Mature" mass to its "Nascent" mass.

* **Nascent Evidence ($E_{nascent}$):** This is tracked natively at the nRNA level now (e.g., unspliced and intronic reads).
* **Mature Evidence ($E_{mature}$):** Spliced exonic reads are tracked at the *transcript* level. To get the mature mass of an nRNA, we simply **sum** the mature evidence of all transcripts that belong to it.

### 2. The Python Implementation

Here is the fully rewritten function. You can drop this directly into `estimator.py`. Notice how much cleaner it is without the `tss_group` logic, and notice how it strictly adheres to our previous breakthrough: **the final EM pseudo-count is a constant $5.0$.**

```python
import numpy as np
import pandas as pd

def compute_nrna_frac_priors(
    n_nrna: int,
    t_to_nrna: np.ndarray,             # [n_t] mapping transcript -> nRNA
    t_mature_counts: np.ndarray,       # [n_t] exonic spliced reads per transcript
    nrna_nascent_counts: np.ndarray,   # [n_nrna] unspliced + intronic reads per nRNA
    nrna_locus_ids: np.ndarray,        # [n_nrna]
    nrna_strands: np.ndarray,          # [n_nrna]
    nrna_global_prio: float = 2.0,
    kappa_locus: float = 10.0,
    kappa_nrna: float = 5.0            # Constant pseudo-counts fed to the EM
) -> tuple[np.ndarray, np.ndarray]:
    """
    Computes Empirical Bayes priors for the nascent fraction (eta) of each nRNA.
    Hierarchy: Global -> Locus-Strand -> nRNA.
    """
    
    # 1. Aggregate mature evidence from transcripts up to their parent nRNA
    nrna_mature_counts = np.zeros(n_nrna, dtype=np.float64)
    np.add.at(nrna_mature_counts, t_to_nrna, t_mature_counts)
    
    nrna_total_counts = nrna_mature_counts + nrna_nascent_counts

    # 2. Global Level (The ultimate fallback)
    global_nasc = np.sum(nrna_nascent_counts)
    global_total = np.sum(nrna_total_counts)
    
    global_eta = 0.5
    if global_total > 0:
        global_eta = global_nasc / global_total

    # 3. Locus-Strand Level Shrinkage
    # Group nRNAs by locus and strand to find the local neighborhood average
    df = pd.DataFrame({
        'locus': nrna_locus_ids,
        'strand': nrna_strands,
        'nasc': nrna_nascent_counts,
        'total': nrna_total_counts
    })
    
    # Sum counts per locus-strand
    grouped = df.groupby(['locus', 'strand'], as_index=False).sum()
    
    # Shrink locus-strand fractions toward the global fraction
    grouped['locus_eta'] = (
        (grouped['nasc'] + kappa_global * global_eta) / 
        (grouped['total'] + kappa_global)
    )
    
    # Map the shrunk locus_eta back to each individual nRNA
    df = df.merge(grouped[['locus', 'strand', 'locus_eta']], on=['locus', 'strand'], how='left')
    locus_eta_arr = df['locus_eta'].values

    # 4. nRNA Level Shrinkage (The final estimate)
    # Shrink individual nRNA fractions toward their locus-strand parent
    nrna_eta_shrunk = (
        (nrna_nascent_counts + kappa_locus * locus_eta_arr) / 
        (nrna_total_counts + kappa_locus)
    )

    # 5. Convert to Beta(alpha, beta) parameters for the EM solver
    # We use a CONSTANT kappa_nrna (e.g., 5.0) so the prior doesn't overpower the EM
    alpha_out = nrna_eta_shrunk * kappa_nrna
    beta_out = (1.0 - nrna_eta_shrunk) * kappa_nrna

    return alpha_out, beta_out

```

### Why this is mathematically beautiful for your pipeline:

1. **No Double Counting:** Because you aggregate the mature transcripts *up* to the nRNA *before* calculating the fractions, you naturally account for all splicing children without double-counting the nascent reads.
2. **Biological Smoothing:** If a gene has a massive spike of intronic coverage on the positive strand, `kappa_locus` will gently pull the prior of a low-coverage nRNA *on that same strand* upwards, giving the EM a highly intelligent starting point.
3. **The `kappa_nrna = 5.0` Shield:** As we discussed in our previous deep dive, feeding `kappa_nrna` as a gentle constant of 5 pseudo-counts ensures that the EM solver is guided by the prior, but is completely free to overwrite it if the iterative read assignment proves the prior wrong.

You don't need to track `tss_groups` anymore. The genomic spans of your decoupled nRNAs have already solved that clustering problem for you perfectly.


========

You caught a major conceptual regression in my previous proposal.

By tossing out the `density + strand + gDNA` hybrid model in favor of raw counts, my proposal fundamentally ignored the physical realities of RNA-seq (e.g., that a 100kb intron will naturally attract 50x more background reads than a 2kb exon, and that gDNA is double-stranded).

Your proposed **Hybrid Plan**—combining the structural elegance of the 3-tier `nRNA -> Locus -> Global` hierarchy with the mathematical rigor of the density/strand model and MoM estimation—is the optimal solution.

Here is my critique and refined implementation plan, resolving your open questions and detailing the exact NumPy data flows required to make this work.

---

## 1. Executive Summary & Critique

The objective is to compute Bayesian Empirical priors for the nascent RNA fraction ($\eta$) to feed into the newly decoupled EM solver.

**The Structural Win:** The old 4-tier hierarchy (`Transcript -> TSS Group -> Locus -> Global`) is collapsed into a 3-tier hierarchy (`nRNA -> Locus-Strand -> Global`). Because an nRNA is defined by its genomic span (start to end), it naturally supersedes and replaces the concept of a "TSS Group".

**The Mathematical Win (The Hybrid Model):** We retain the rigorous `_compute_hybrid_nrna_frac_vec()` function. This preserves:

1. **Density Normalization:** Adjusting for the massive length differences between introns and exons.
2. **gDNA Subtraction:** Removing the global background gDNA density from the intronic signal.
3. **Strand Correction:** Utilizing the inverse-variance $(2s-1)^2$ weighting to seamlessly scale between unstranded and stranded libraries.

## 2. Resolution of Open Design Questions

### Q1: How to compute nRNA-level exonic length?

*Your options: (a) mean of member transcripts, (b) max, (c) actual union geometry.*
**Recommendation:**
We should compute this as a weighted average of the effective lengths of the nRNA's "child" transcripts. A nascent RNA can give rise to multiple mature RNAs (splicing isoforms). Each of the transcripts will have estimated counts and effective length. The effective length of the group of mature RNAs can be computied as weighted average of the nascent RNA's mature RNAs using the current estimated abundance of each mature RNA. This way, if one of the transcripts dominates (ex. >99% of counts, 2kb effective length), then that transcript will dominate the exonic effective length calculation.

### Q2: Keep MoM or go constant $\kappa$?

**Recommendation:** Keep MoM to estimate `kappa_global` and `kappa_locus` because dataset variance fluctuates wildly, and MoM calculates exactly how much to trust the local biological baseline.
However, for the *final* output fed to the EM solver, use a constant `kappa_nrna`.
*Why?* The MoM establishes the *prior mean* ($\eta_{shrunk}$). Multiplying this mean by a gentle, constant ensures the EM is always guided by the baseline but is never paralyzed by an overwhelmingly heavy prior if the local EM evidence contradicts it.

### Q3: Keep `_compute_hybrid_nrna_frac_vec()`?

**Recommendation:** Yes, keep it intact. It is beautifully factored. We just change the inputs from Transcript-level arrays to nRNA-level arrays.

---

## 3. Refined Implementation Plan

### Phase 1: nRNA-Level Data Aggregation

To feed `_compute_hybrid_nrna_frac_vec()`, we must bridge the transcript-level exonic data to the nRNA-level intronic data.

**1. Exonic Counts (Aggregated):**

```python
nrna_exonic_sense = np.zeros(n_nrna, dtype=np.float64)
np.add.at(nrna_exonic_sense, t_to_nrna, transcript_exonic_sense)

nrna_exonic_anti = np.zeros(n_nrna, dtype=np.float64)
np.add.at(nrna_exonic_anti, t_to_nrna, transcript_exonic_antisense)

```

**2. Intronic Counts (Native):**
In the new decoupled architecture, `nrna_intronic_sense` and `nrna_intronic_antisense` are already natively collected at the nRNA level during the BAM scan. No aggregation needed.

**3. Lengths:**

Weighted average of transcript exonic lengths. The extreme case is if only one of the nrna group of transcripts is expressed, and the rest of the transcripts are ZERO. We don't want the mean of effective lengths in this case -> we just want the effective length of the dominant transcript isoform. So a count-weighted average should give us an appropriate nrna 'exonic' effective length estimate. Consider the previous case where we each transcript has it's own nRNA shadow. In that case, we were computing nRNA exonic length based on the transcript effective length. Now we are computing nRNA exonic length based on a group of transcripts. Theoretically the count-weighted average should be the mathematically equivalent, is that right? 

# Intronic length is simply the span minus the mean exonic length
nrna_intronic_length = np.maximum(1.0, nrna_spans - nrna_exonic_length)


**4. Locus IDs:**

Nascent RNAs share the same locus as their children transcripts

```python
# Since all transcripts in an nRNA share a locus, we can just grab the first one
nrna_locus_ids = transcript_locus_ids[nrna_to_t_offsets[:-1]]

```

### Phase 2: The 3-Tier Shrinkage Pipeline

Rewrite `compute_nrna_frac_priors()` to execute the 3-tier hierarchy.

1. **Base Estimates:** Call `_compute_hybrid_nrna_frac_vec()` on the arrays generated in Phase 1. This yields `nrna_raw_eta` and `nrna_variance`.
2. **Locus-Strand Estimates:** Group the arrays from Phase 1 by `(nrna_locus_ids, nrna_strands)`. Sum the counts/lengths and call `_compute_hybrid_nrna_frac_vec()` on the aggregated locus-strand vectors to get `locus_raw_eta`.
3. **Global Estimate:** Calculate the global evidence-weighted mean of `locus_raw_eta`.
4. **MoM Shrinkage:** * Shrink `locus_raw_eta` toward `global_eta` using MoM `kappa_global`.
* Shrink `nrna_raw_eta` toward `locus_shrunk_eta` using MoM `kappa_locus`.
* *Result:* `nrna_eta_shrunk`.


5. **Final EM Prior:**
```python
estimator.nrna_frac_alpha = nrna_eta_shrunk * kappa_nrna_constant  # e.g., 5.0
estimator.nrna_frac_beta = (1.0 - nrna_eta_shrunk) * kappa_nrna_constant

```



### Phase 3: Config and Dead Code Eradication

* **Delete:** `compute_tss_groups()` in `index.py`, `t_to_tss_group`, and `_aggregate_nrna_frac_by_group()`.
* **Config Updates:** * Remove `kappa_tss`, `tss_window`, `mom_min_evidence_tss`.
* Keep `nrna_frac_kappa_global`, `nrna_frac_kappa_locus`.
* Add `nrna_frac_kappa_nrna: float = 5.0` (The constant prior strength fed to the EM).



### Phase 4: Callers & Tests

* **`pipeline.py` & `profiler.py`:** Update to call the new `compute_nrna_frac_priors` directly with the nRNA-level locus/strand arrays.
* **`test_estimator.py`:** * Delete TSS-group specific tests.
* Add assertions ensuring that `alpha + beta == kappa_nrna_constant` (verifying the final EM constraint).
* Update golden datasets (expect minor shifts in $\eta$ due to the shift from transcript-averaging to native nRNA computation).



---

### Final Thoughts on the Plan

Your pushback saved the mathematical integrity of the nascent estimation. By structuring it exactly as outlined above, you preserve the precise handling of fragment lengths, gDNA density, and protocol strandedness, while still cashing in on the massive structural simplifications of the decoupled nRNA EM solver.