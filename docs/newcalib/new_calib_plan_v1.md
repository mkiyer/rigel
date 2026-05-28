Your proposal to transition from a discrete, classification-based hybrid capture model to a continuous, Empirical Bayes (EB)-shrunk **exposure factor framework** is theoretically sound, elegant, and addresses several core pain points in the current pipeline architecture.

Below is a detailed, structured audit of your proposed redesign, covering theoretical considerations, mathematical validation, code-level impacts, and critical bug fixes for downstream execution.

---

### 1. Theoretical Audit: The Continuous Exposure Concept

**Verdict: Highly Recommended. This is the correct statistical direction.**

Hybrid capture enrichment is inherently continuous, not binary. Factors such as bait tiling density, GC content, local secondary structures, and hybridization kinetics create a wide, highly skewed spectrum of capture efficiency (ranging from complete dropouts to $>1000\times$ amplification).

* **The Problem with the 4-State Model:** Forcing this continuous spectrum into discrete latent strata (`STATE_UNEXPRESSED_CAPTURE` vs. `STATE_UNEXPRESSED_OFFTARGET`) introduces arbitrary classification thresholds. This creates severe instability in low-coverage regions, makes the model highly vulnerable to extreme class imbalances (e.g., small panels targeting a handful of genes), and forces the system to try to globally classify whether a library is "capture-on" or "capture-off."
* **The Continuous Exposure Advantage:** Modeling capture depth as a localized continuous parameter—an *exposure factor* ($\omega_r$) representing relative sampling efficiency—unifies both capture and non-capture data under a single framework. For non-capture datasets, the exposure factors naturally stay tightly bounded around $1.0$ (uniform background). For capture datasets, the distribution expands organically to accommodate extreme localized scaling without requiring hard thresholds or separate pipelines.

---

### 2. Deep Interrogation of Current Code & gDNA Mass Estimation

In the current v6 architecture, the calibration iteration (`run_calibration_iteration`) performs a 4-state EM loop using a vectorized log Bayes factor tensor.

* **Reliability of Current Estimators:** The pipeline already contains highly capable foundational blocks for estimating raw gDNA mass:
1. *Stranded Data:* `deconvolve_compartments_by_strand` cleanly isolates antisense RNA expression from unstranded genomic contamination to yield a reliable `RegionGdnaChannelEstimate`.
2. *Unstranded/Intergenic Data:* `predict_contained_gdna_from_excess` and `run_boundary_sweep` utilize boundary-crossing fragments to impute contained gDNA mass where structural information is sparse.


* **What Must Change:** If you eliminate the discrete capture states, the 4-way classification tensor disappears. The EM loop collapses back to a simple, highly robust 2-state model: **Expressed vs. Unexpressed**. You will use the posterior probability of the region being *unexpressed* ($P(\text{expressed}) \approx 0$) to isolate clean, unexpressed genomic mass across the library. This unexpressed mass is then compared directly against the global background rate ($\rho_{\text{off}}$) to calculate the raw, unshrunk exposure ratios.

---

### 3. Native Fractional Accumulator Code & Discrete Count Instrumenting

Your intuition regarding fractional mass vs. discrete observations is exactly correct, though your overlap mathematics contains a minor offset.

* **Mathematical Correction:** For a genomic region with boundaries $(a, b)$ and width $W = b - a$, and a fragment of physical length $L$, the range of potential fragment start positions that can overlap this region spans from $(a - L + 1)$ to $(b - 1)$. The total number of valid start positions is:

$$\text{Opportunity} = (b - 1) - (a - L + 1) + 1 = (b - a) + L - 1 = W + L - 1$$



For a tiny $1\text{ bp}$ region ($W=1$), there are exactly $L$ starting positions where a fragment can overlap it. Therefore, even a $1\text{ bp}$ region can have thousands of physical fragments contributing fractional scores to it if sequencing depth is high.
* **The Shrinkage Problem:** Currently, the native accumulator tracking only stores accumulated fractional mass. This is insufficient for Empirical Bayes because fractional mass obscures the **Effective Sample Size (ESS)**. A regional gDNA mass of $0.2$ derived from a single fragment with $20\%$ overlap has high variance; a mass of $0.2$ accumulated from 200 separate fragments with $0.1\%$ overlap each has extremely low variance.
* **The Fix:** You must instrument the native C++ accumulator to track an integer count of unique physical fragment observations ($N_r$) touching each region. During the single-pass BAM scan, every time a fragment intersects region $r$, increment $N_r$ by $1$. This provides a discrete count that acts as the statistical weight for your shrinkage algorithm.

---

### 4. Moving 'Boundary' Code to Native C++

**Verdict: Highly recommended to execute now.**

Currently, `BoundaryTable` acts as transitional Python-side storage, splitting up data after the native accumulator pass. Since you already need to modify the native C++ structures to capture discrete fragment counts (Item 3), this is the optimal time to move the region/boundary structural definitions into native storage. Assigning fragments directly to either `ContainedRegion` or `Boundary` slots during the single-pass C++ BAM scan will reduce peak memory footprint and eliminate the overhead of Python-side dictionary and tensor reorganization.

---

### 5. Downstream Execution: The Denominator Normalization Bug

You have identified a **critical, severe bug** in how exposure factors are handled downstream.

* **The Current Broken Behavior:** In `gdna_prior_current_construction.md`, the pipeline computes `gdna_eff_len[locus] = FL-marginal effective length * mean regional A_r`. In the native solver, this effective length enters the denominator of the read likelihood ($P(\text{read} \mid t) \propto 1/L_{\text{eff}}$).
* **Why this destroys performance:** If a region is highly enriched by capture (e.g., $A_r = 1000$), multiplying the effective length by $1000$ makes the component look structurally massive to the EM solver. This drastically *reduces* the likelihood of individual fragments mapping to the gDNA component, penalizing it in the exact regions where capture is highly efficient. The gDNA mass is artificially suppressed and "siphons away" into the competing RNA components.
* **The Correction:** Exposure ($\omega_r$) scales *sampling visibility*. Enriched regions produce more reads for the same physical mass. Therefore, in the generative model, a high exposure factor should *increase* read likelihood. The effective length parameter passed downstream must be **divided** by the exposure factor, not multiplied:

$$L_{\text{eff\_adjusted}} = \frac{L_{\text{eff}}}{\omega_r}$$



This shrinks the effective length penalty for on-target loci, making the gDNA component highly competitive against RNA in captured regions, matching physical reality.

---

### 6. Big Picture: Designing the EB Shrinkage Framework

To implement your vision of a continuous exposure factor that shrinks gently to the global center based on data sparsity, use a **Gamma-Poisson / Negative Binomial hierarchical framework**.

#### The Mathematical Formulation

Let:

* $\hat{M}_r$ = Estimated continuous gDNA mass in region $r$ (derived from strand deconvolution or boundary imputation).
* $L_r$ = Base effective length of region $r$.
* $\rho_0$ = Global background baseline density (estimated from high-confidence unexpressed, uncaptured seed regions).
* $\lambda_r = \rho_0 \cdot L_r$ = The expected uniform background count for the region.

We place a Gamma prior over the true continuous exposure factor $\omega_r$:


$$\omega_r \sim \text{Gamma}(\nu, \nu)$$


Setting both shape and rate to $\nu$ enforces a prior mean of $1.0$ (representing uniform, non-capture global sequencing opportunity). The parameter $\nu$ acts as your **prior effective sample size (ESS)**, determining the gravitational pull toward the global center.

Updating this prior with the observed regional evidence yields a conjugate posterior distribution for the exposure factor:


$$\omega_r \mid \hat{M}_r \sim \text{Gamma}(\nu + \hat{M}_r, \nu + \lambda_r)$$


The EB-shrunk exposure factor is the posterior mean:


$$\bar{\omega}_r = \frac{\nu + \hat{M}_r}{\nu + \lambda_r}$$

#### Dynamic Variance-Aware Shrinkage (The Solution to Brittleness)

Instead of hardcoding a global $\nu$ (which is brittle), **learn $\nu$ from the sample variance of the data**.

1. Isolate all regions that are highly confident unexpressed regions ($P(\text{expressed}) \approx 0$).
2. Compute the raw, unshrunk ratios $R_r = \hat{M}_r / \lambda_r$ for these regions.
3. Calculate the sample variance $\sigma^2$ of these ratios across the genome.
4. Set $\nu = \frac{1}{\sigma^2}$.

#### How This Handles Your Scenarios Automatically

* **Non-Capture Datasets:** The data is well-behaved and uniform; the variance $\sigma^2$ of the ratios will be extremely small. Thus, $\nu$ will be very large (e.g., $\nu = 100$). The prior exerts a powerful gravitational pull, forcing all regional exposures to remain tightly bound around $\bar{\omega}_r \approx 1.0$.
* **Capture Datasets:** The enrichment profile is highly skewed, generating a massive sample variance $\sigma^2$. Consequently, $\nu$ will collapse to a tiny value (e.g., $\nu = 0.05$). The prior gravity vanishes, allowing sparse or heavily enriched regions to dictate their own exposures completely unhindered:

Let's look at your examples under a capture setting where variance is high and $\nu = 0.1$, given a global background expectation of $\lambda_r = \rho_0 \cdot L_r$:

* **10bp region, 2 fragments ($\hat{M}_r = 2$, $\lambda_r = 0.01$):**

$$\bar{\omega}_r = \frac{0.1 + 2}{0.1 + 0.01} = \frac{2.1}{0.11} \approx 19.09 \quad (\text{Maintains strong local enrichment profile despite small size})$$


* **100bp region, 20 fragments ($\hat{M}_r = 20$, $\lambda_r = 0.1$):**

$$\bar{\omega}_r = \frac{0.1 + 20}{0.1 + 0.1} = \frac{20.1}{0.2} \approx 100.5 \quad (\text{Allows data to take over as count size scales up})$$


* **100bp region, 200 fragments ($\hat{M}_r = 200$, $\lambda_r = 0.1$):**

$$\bar{\omega}_r = \frac{0.1 + 200}{0.1 + 0.1} = \frac{200.1}{0.2} \approx 1000.5 \quad (\text{Fully captures the massive } >1000\times \text{ capture spikes without clipping})$$



### Summary Action Plan

1. **Collapse the latent states** from 4 to 2 (Expressed vs. Unexpressed) in `latent_states.py` and `calibration_iteration.py`.
2. **Instrument the native C++ BAM scanner** to track integer fragment counts $N_r$ alongside fractional mass accumulation, and port the transitional `BoundaryTable` logic into native memory space.
3. **Implement the Variance-Informed Gamma EB Shrinkage** using the unexpressed regions to dynamically calculate the prior strength $\nu$.
4. **Deploy a downstream PR** to correct the denominator normalization bug by changing $L_{\text{eff}} \cdot \omega_r$ to $L_{\text{eff}} / \omega_r$.