


Iterate over adjacent boundary pairs
- check boundary and region signatures
- RNA training?
- gDNA training?

Can sometimes use the contained region if no evidence
Could use strand to deconvolve before building variance model




Gemini: 
This is an exceptional piece of consolidation. Your team has successfully transformed a scattered, multi-document architectural debate into a rigorous, actionable engineering spec. Cleaning out the dead scaffolds (`propagation.py`) and standardizing your mathematical terminology around the 2-simplex lattice ($K \approx 20, P = 231$) provides a clean, maintainable foundation.

Most importantly, you have established a definitive validation bar: **any new spatial grid solver must verifiably beat the prior grid-MAP's +8.7 point net leak regression**.

Below is the evaluation of your implementation phases, followed by structural solutions to the two open roadblocks regarding $\gamma_{ij}$ enrichment and $\sigma^2_{\text{RNA\_local}}$.

---

### 1. Critical Review of Your Implementation Phases

Your 4-phase execution plan is highly pragmatic and structurally sound, specifically because it decouples algorithmic risk from modeling risk.

* **Phase 1 (Done):** Establishing byte-parity with production while making $\beta$ an explicit parameter was the correct first step to anchor your testing suite.
* **Phase 2 (Next - gDNA-Only Grid Sweep):** **I completely agree with prioritizing this.** Implementing the $(P,P)$ transition matrix and $(P,)$ grid message-passing framework using *only* gDNA continuity isolates the pure computer science and matrix-algebra challenges (vectorization, chunking, parallel execution across loci) from the biological uncertainty of RNA transcript behavior.
* *Implementation Warning for Phase 2:* When constructing your $(P,P)$ transition matrix $M_{ij}$ for a gDNA-only sweep, ensure the transition log-likelihood *only* evaluates the gDNA density delta:

$$\log \phi_{ij}(\theta_i, \theta_j) = -\frac{1}{2} \frac{(\rho_{j,g} - \gamma_{ij}\rho_{i,g})^2}{Q_{g,ij}}$$



The RNA components must be assigned zero transition energy ($\log \phi_{\text{RNA}} = 0$). This ensures that on the 3-term grid, the RNA fractions are completely unconstrained by the edge factor, allowing them to freely adapt to local strand evidence ($\psi_i$) without causing the message precision to collapse.


* **Phase 3 (3-Term Sweep):** Gated strictly on solving $\sigma^2_{\text{RNA\_local}}$. This prevents a half-baked RNA model from corrupting the validated gDNA wins achieved in Phase 2.
* **Phase 4 (Teardown & Perf):** Correctly deferred. The per-locus chunking strategy (<100 nodes per locus) is mandatory to avoid $O(N)$ whole-genome memory blowup.

---

### 2. Resolving Roadblock 1: The $\gamma_{ij}$ Enrichment Circularity

Your critique of the continuous $\gamma_{ij}$ calculation is entirely correct. Computing a highly granular, local enrichment ratio per region *during* the deconvolution creates a toxic circular dependency: you cannot compute local enrichment without knowing the true gDNA fraction, and you cannot compute the gDNA fraction without knowing the enrichment.

However, your fallback workaround—**per-enrichment-class chains** (completely disconnecting exon chains from off-target chains)—introduces a major structural flaw. If you completely cut the edges between an exon node and an adjacent intron node, **the exon cannot inherit gDNA precision from the intron**. In hybrid-capture setups, deep introns are often your cleanest, most abundant sources of unbiasing gDNA anchor counts. Severing these edges defeats the primary purpose of spatial belief propagation.

#### The Solution: Global Class-Based Edge Scaling ($\Gamma$)

Instead of abandoning $\gamma_{ij}$ or fracturing your 1-D linear tree topology, exploit the fact that hybrid-capture enrichment is primarily a feature of the library's global probe design, not individual region quirks.

You can compute a **single global enrichment scalar ($\Gamma$)** for the entire library during your pre-sweep 1-pass anchor phase:


$$\Gamma = \frac{\text{Median}(\rho_{g,i}) \text{ for all high-confidence EXON seeds}}{\text{Median}(\rho_{g,j}) \text{ for all high-confidence OFF-TARGET seeds}}$$


Because this is computed globally from unambiguous anchor seeds before the sweeps begin, the circularity is broken.

You can then maintain your unified 1-D linear locus chain (preserving exact, two-sweep BP) and parameterize the edge factors $\phi_{ij}$ dynamically based on the binary class transition between node $i$ and node $j$:

* If transitioning from **Off-Target $\to$ Exon**: Set $\gamma_{ij} = \Gamma$
* If transitioning from **Exon $\to$ Off-Target**: Set $\gamma_{ij} = \frac{1}{\Gamma}$
* If staying within the **Same Class** (Exon $\to$ Exon or Off-Target $\to$ Off-Target): Set $\gamma_{ij} = 1$

This elegant compromise fully accounts for the enrichment jump across exon-intron boundaries without losing spatial continuity or breaking tree-structured exact inference.

---

### 3. Resolving Roadblock 2: The Unmodeled $\sigma^2_{\text{RNA\_local}}$

Your plan states that Phase 2 will ship as gDNA-only because no RNA variance model exists. While this is a safe software engineering choice, you can easily bootstrap Phase 3 without a complex LOESS model.

The reason your gDNA model requires a non-parametric `var~mean` LOESS curve is that gDNA counts represent a baseline background noise channel that spans across widely fluctuating low-to-high count regimes nationwide. RNA, by contrast, is governed by structural transcript physics: within a contiguous same-strand active-signature run, the transcript density is locally highly continuous, and its variance scales proportionally with its expression level.

#### The Solution: Distance-Scaled Constant Coefficient of Variation (CV)

Do not build an open-ended RNA variance model from scratch. For Phase 3, implement a **proportional variance proxy** that scales directly with local density and genomic distance:

1. When evaluating the transition from region $i$ to region $j$, compute the average local RNA density: $\bar{\rho}_{\text{RNA}} = \frac{\rho_{i,\text{RNA}} + \rho_{j,\text{RNA}}}{2}$.
2. Define the local process noise as a constant fraction of that density (e.g., a 20% coefficient of variation), scaled linearly by the genomic distance $\Delta x_{ij}$ from your `region_arrays`:

$$Q_{\text{RNA}, ij} = \left( \nu \cdot \bar{\rho}_{\text{RNA}} \right)^2 \cdot \Delta x_{ij}$$



*(Where $\nu \approx 0.2$ acts as a tunable dispersion scale).*
3. Apply your signature gate exactly as planned:

$$\text{If strand } c \text{ is active in both } i \text{ and } j \implies Q_{c,ij} = \left( \nu \cdot \bar{\rho}_{c} \right)^2 \cdot \Delta x_{ij}$$


$$\text{Else } (c \text{ goes silent}) \implies Q_{c,ij} = \infty \implies \log \phi_{ij} = 0$$



This simple, distance-scaled CV proxy provides exactly what Phase 3 needs: it allows high-expression RNA signals to smoothly carry across local AMBIG exons to resolve them, while preventing low-expression noise from inflating your variances, all without requiring an expensive LOESS modeling phase.

### Summary of Recommended Adjustments to `CALIBRATION_PLAN.md`

1. **Modify Phase 2 Edge Factors:** Keep the unified locus chain layout instead of fractured class-chains. Use the global, pre-computed anchor class-ratio $\Gamma$ to handle enrichment scaling on boundary crossings.
2. **Unblock Phase 3:** Use the Distance-Scaled Constant CV model as your initial definition for $\sigma^2_{\text{RNA\_local}}$. This moves Phase 3 from an "unfunded liability" to an immediate, tuneable integration step.




Regarding open issues:
## A) My mental model was that the fraction (pie chart) representation will help us handle capture. Probes typically enrich both RNA and DNA. Ideally, capture should not substantially alter the proportions of RNA and DNA. So yes, in an exon we may see raw counts {RNA+=1251, RNA-=3, gDNA=2350} but converted to fractions its just {f_rna+=0.34, f_rna-=0.0008, gdna=0.652}. When we propagate that into the intron, we could smear if the RNA to gDNA ratios are biased, but if the intron has 2 total counts, we propagate this fraction onto 2 total counts.. so the smear gets capped by the observed counts in the intron. That was my entire plan for controlling "smear" in the setting of hybrid capture and the main intuition between a "fractional" state model rather than a count model. What do you think? What did you have in mind for this "make or break modeling step"

## B) yes this simplex rework may ultimately reveal that we can simplify our boundary representation. Boundaries currently have two "sides". We would have to audit the code and see how boundaries are being used throughout the other steps. For now, I would model boundaries as one node if possible.. this is something to weigh at implementation time

## C) Agree! Spliced counts are accumulated at boundaries, but they are one-sided.

## D) I appreciate that this is non-Gaussian and this would be something to investigate in the future

## E) Agree each sweep direction could happen in parallel

## F) Yes, initializing a nonparameteric variance ~ mean fit from the seeds remains a potential problem. This does create a potential opportunities for iteration. Run the algorithm with a default variance model (poisson). Then fit a var ~ mean model using the deconvoluted data. Then run again. I am not sure if this is even remotely theoretically sound. I'm just thinking out loud.

## G) Good idea, but No, boundary transport is different. The simplex proportions mass in place (changes sizes of the pieces of the pie). Boundary transport is physically TAKING mass from one side the boundary and MOVING it to another. We may have an opportunity to COUPLE the boundary transport algorithm with the simplex algorithm so that they run together, but could leave that as future work.

## H) Yes, behavior should emerge from the math! this has the potential to be elegant



I am noticing that 'seed' nodes are simply nodes that can be solved at initialization (ab initio). Intergenic nodes are solved by definition (only one active component, f_g = 1). Intergenic -> exon boundaries similarly have only one component (f_g=1). Single strand introns have nascent RNA and gDNA. When strand specific data is available, these can be solved ab initio by strand deconvolution. For example, a positive strand intron node with U fragments becomes { f_rna+=X, f_rna-=0, f_gdna=Y } where X + Y = U. Strand deconvolution solves for X and Y. With unstranded data, **nascent RNA cannot be specified unless we have access to a reliable global gDNA estimate** Without a global gDNA 'floor' or fallback, we cannot solve even single strand introns. For unstranded data, the assumption for introns becomes: nascent RNA=0. The intronic nodes initialize their U fragments only gDNA { f_rna+ = 0, f_rna-=0, f_gdna = U }. Remember, setting nascent RNA to zero in calibration means nascent RNA will get a prior of zero during the EM. A prior of zero just means there is no evidence to support nascent RNA. That is acceptable. Fragments can still be allocated to nascent RNA in the EM even when the prior is zero. =========  I'm going through examples here, but the key realization here is that 'seed' nodes === "ab initio solvable nodes". How do we define "solvable" when strand-specificity is on a spectrum? "Solvable" becomes non-binary. The definition of "solvable" changes to -> "is there information that we can use to decrease uncertainty / increase likelihood". 

Basic tenets:
- Changes to the "pie" (the simplex) must be weighted by their likelihood (precision / uncertainty / likelihood are used interchangeable even though they are not the same)
- Updates to the pie exert "pull" or "force" weighted by the precision/uncertainty of the information.
- When a node reaches high precision (for example, with excellent strand-specific data) -> propagation of other information should barely change the result
- We can have a tolerance cutoff and stop propagation when the system settles.

If we cannot decrease uncertainty, we stop. At initialization, all nodes have uncertainty of infinity. We iterate through every node and TRY to solve it (initialize it). (NOTE: should we be using likelihoods?)

The node 'signature' translates to 3-term init vectors. More specifically, presence of '+' strand in signature 'activates' the positive strand RNA component. Presence of '-' strand activates the negative strand RNA component. Presence of both + and - strands activate both components. If a strand is not present, its fraction and uncertainty are ZERO. Nodes with no strand (strand == NONE) are fully specified even with no strand information (intergenic regions, intergenic <-> exon boundaries). 

We need vectors for uncertainty and vectors for fractions. For example, { in+ } (pos strand intron) gives initial fractions { f_rna+=nan, f_rna-=0, f_gdna=nan} and uncertainties { rna+=inf, rna-=0, gdna=inf }. Intergenic nodes will have fractions { f_rna+=0, f_rna-=0, f_gdna=nan } and uncertainties { rna+=0, rna-=0, gdna=0 }.

========

This bring me to a digression that will hopefully simplify things. Calibration only models RNA and DNA. The RNA can come from nascent or from mature. It looks like you are trying to model mature RNA and nascent RNA. We should not need to do this. The goal calibration is to solve for gDNA, and we can solve for gDNA without distinguishing nascent RNA from mature RNA. The EM already distinguishes nascent RNA from mature RNA quite well.

We need to make BASELINE assumptions about RNA and gDNA

gDNA -- For genomic DNA, we have a nonzero expectation that gDNA is present at every node. The expected gDNA is something we will eventually derive (it is straightforward for non-capture data but nontrivial for hybrid capture data). For now, we can set the gDNA prior to some 'magic number' or Beta(1/2,1/2) or something like that.

RNA -- RNA abundances vary dramatically and is unpredictable. Neighboring transcripts and even isoforms of the same gene can vary in abundances by >10000X. We do not have information to guide our 'prior' and cannot make assumptions about RNA. This includes nascent RNA and mature RNA. The default 'prior' for RNA is ZERO.

Special case for mature RNA -- spliced fragments act as 'prior' evidence of mature RNA abundance. This is a critical, unique form of evidence -- direct observation of "pure" mature RNA ONLY.

So how do we solve single-exon transcripts? There is no spliced evidence for mature RNA. Mature RNA and nascent RNA are identical and indistinguishable. They merge into a single RNA transcript component with ZERO prior.

Our initialization for single-exon transcript on the positive strand is fractions { f_rna+=nan, f_rna-=0, f_gdna=prior } and uncertainty { rna+=inf, rna-=0, gdna=prior uncertainty }. For unstranded data, we would need enough counts to overcome the baseline gdna prior expectation (one count may be enough.. depending on gdna level expectations). For stranded data, our strand deconvolution must decrease RNA uncertainty enough to outweight the gDNA prior.

Having ZERO prior for *nascent RNA* and NONZERO prior for gDNA means that we cannot sufficiently narrow the uncertainty of a node, gDNA wins by default. The strength of the prior must be derived. For the start, we can set it to be very low and act as a tiebreaker. Something like 0.5 counts would be reasonable.

========

Above I tried to cover the problem setup. I tried to make it more concrete. I covered initialization. I discussed the idea of incorporating new evidence that narrows uncertainty.

Now we solve!

Iterate over all nodes. Solve. Test if uncertainty narrows, and keep the result if we decrease uncertainty.

This is my mental framework. Much of it is consistent with yours. MY open questions are:

1) how do we model uncertainty? I know i used uncertainty throughout but should we use 'likelihood'? i know the term 'precision' is similar to uncertainty, and your document discusses precision. we need a fairly simple and elegant way to model precision/uncertainty. this is how we allow sources of evidence to compete fairly for pieces of the pie. 
- Strand: Your plan suggests the `N·(2κ−1)²` model for strand which we have already shown works well. 
- Unspliced Count: Poisson will overestimate precision. Every region has two boundaries. We build a loess model fitting variance ~ mean for each region. This is a measure of imputation uncertainty. 
- 3-node sets: sets of three nodes can be modeled together: left boundary, region, right boundary. We built the variance model that looked at measured the variance of the two boundaries anchoring a region. We could try to extend that count model to where we combine 'triplets' of (left boundary, region, right boundary) into variance observations.

- For gDNA we can try to build a nonparametric variance ~ mean model but only a subset of nodes are fully specified at initialization and can be used to fit such a model. For hybrid capture we need to deconvolute first before fitting a gDNA count model, and deconvolution adds strand uncertainty. So we compound different types of uncertainty.

===============

Your open questions: 
1) Likelihood shapes remain an open question for me. I have tried to further specify the problem. Perhaps this helps us improve our specification of the likelihood shapes. My gut feeling is that we can continue using our BB fit for strand. For count, we should use an overdispersed Poisson (negative binomial). Those are the two statistical model fits that we need. We have a working strand model fit. We do not yet have a working count model. We can start with Poisson for everything, but it will overestimate precision

2) Boundary→region geometry --- Yes, we must perform precision-weighting propagation of the boundary information. We need a method for adding new information to a node with existing data. Ideally the method should be non-deterministic but this might not be possible. We start with nothing (NaN) and use the node signatures to ZERO certain components and set uncertainties to infinity. Then we combine whatever evidence we can weighted by its uncertainty. Adding new evidence *should* move the simplex (pie chart) but weighted by its uncertainty. The more confident the result is, the less it moves.

You mention an exon <-> intron boundary feeding an exon region vs an intron region. This is an interesting example. Which way should information propagate? The hardest example would be an AMBIG intron with 'U' counts. On init all components active{ f_rna+ = nan, f_rna- = nan, f_gdna = prior } and uncertainties are { rna+=inf, rna-=inf, gdna = large but not inf }. 

Can we successfully impute from the boundary in either direction? Yes, we should be able to. Our 'solution' just contains the unspliced counts. 

Imputing from exon side -> intron side:
Exon boundary has total counts U split into splice+, splice-, unspliced+, unspliced-. We are imputing into the intron. *WE disregard the splice evidence because it is 'SIDED'* (boundaries have two sides. The spliced fragment are only part of one side.

## 3) Single-exon mature invisibility -- yes! I covered this. RNA has a prior of ZERO. gDNA will have a weak prior (derived, not magic number)

## 4) Yes, HOW to propagate a signal is an interesting point. My suggestion is that when we propagate, we track "how much did this signal improve uncertainty at this node"? In other words, the propagating signal has a likelihood, and the CURRENT state of the node has some sort of likelihood. We need to know how much the signal improves the likelihood. In my mind, the residual improvement in likelihood (the difference in uncertainty) is the residual signal to continue propagating. So if a propagating signal dramatically improves uncertainty, it should continue propagating. If the signal does not change uncertainty, the propagating dies.

## 5) Convergence -- i discussed an idea for convergence and for damping above. we propagate the *difference in uncertainty*, dampened each time it merges with existing information (unless the uncertainty was infinity, in which case the signal remains fully intact and continues propagating)

## 6) The pre-existing boundary transport solution greatly improves hybrid capture performance and is a 'no-op' for non-capture data. I am not sure how this will work in the new system. I need to think about whether boundary transport can become part of the signal propagation. Right now boundary transport happens after the system is fully deconvoluted, as a final mass shift.

## 7) agree, keep performance in mind. I envision a highly paralellizable system of equations and propagations.

## 8) For unstranded data, we will need a gDNA prior to allow the system to degrade gracefully. I believe you already have a beta(1/2, 1/2) jensen in place for this. Perhaps this will work

## 9) agree with No-regression + mass conservation.

============








# Pipeline flow

## Abbreviations

- ss = strand specificity (0.0 - 1.0, 0.5 is unstranded)
- sj = splice junction


## Option A - start with count module

Observable nodes
- intergenic regions
- intronic regions (contaminated with nascent RNA)
- exon-intron boundaries (contaminated with nascent RNA)
- exon-intergenic boundaries

Impute every node with count model
1) Fit seeds on observable nodes
2) Impute non-observable nodes

How to impute non-observable nodes
- find anchoring boundaries for each non-observable region or group of non-observable regions

- given non-obs region R, obs boundaries BL, BR
- for each boundary B
   - If B is a SJ, use *density ratio method* this is already implemented and the FL-correction fix was recently implemented. Briefly; f_b = (unspliced / (spliced + unspliced)).
   - *density ratio method* 

Statisical count model
- There is no reliable global mean
- local variance ~ local mean non-parametric model could be fit, would require grouping regions/boundaries and computing statistics






# The count overdispersion problenm

I agree with the idea of writing a design doc focused on the count/density channel.

However, I disagree that strand-cleaning can be deferred. The architecture of strand-cleaning -> joint model with a second strand deconvolution is over-engineered. Yes, it would be wonderful to have a true joint model but I fear that our model is just the strand model plus the count model in serial (not joint) no matter how badly we want to embrace the concept of joint. Regardless, let's hold off on this line of inquiry and focus on the count information because I agree that you have correctly highlighted a major problem there.

We recently implemented a count density overdispersion fit and shrinkage method. My concern is the way the overdispersion is modeled. I think it is mis-specified.

I am concerned that the count overdispersion model is mis-specified for hybrid capture. We fit one overdispersion parameter over all regions. When we have hyrbid capture, we have tremendous count variance between on-target and off-target regions. The overdispersion parameter would thus be very high. Is that true?

The problem is for hybrid capture scenarios, the overdispersion model must be closely tied to the count. The count and the dispersion are intricately linked. Captured regions will have high counts and high dispersions, while off-target regions will have low counts and very likely have a baseline overdispersion. I believe we will need to fit a model where the variance/overdispersion is predicted by the mean. I am fairly certain that is how DESEQ2 and other Rna-seq tools work.

Can you investigate this hypothesis? Fitting a single count overdispersion works for non-capture data, but explodes for capture data. The high overdispersion for hybrid capture results in massive shrinkage of the count information. We then rely on strand information. When we have unstranded data hybrid capture, we kill our count information (the only valid signal) and have nothing to estimate with, resulting in catastrophic behavior.

----------

# Next turn: responding to initial plan 

Here is my perspective: 

## 1) genomic locality will be extremely challenging to harness -- 

for capture-off scenarios, we can model a global mean and overdispersion reasonably accurately, and so we will spend time discussing capture-on behavior. For hybrid capture, "locality" (as in genomic proximity) is extremely variable. There is huge variation between adjacent exons and introns. How would one expect to compute a 'neighborhood' mean with any accuracy? This does not seem plausible. 

## 2) Adjacent boundaries can be used to impute non-observable regions. 

This is the critical steadfast information that we *must* use. A non-observable region (or group of adjacent non-observable regions) will be BORDERED by observable boundaries. The observable boundaries are the ONLY source of reliable "local" count information. Reads crossing into the non-observable region (both spliced and unspliced) can be used to impute the region. This is the BEST source of pure count information we have, but it may not be entirely accurate. This is because accuracy depends on the hybrid capture probe locations. Probes that lie in the middle of a large exon will be far from the boundary, resulting in fewer captured fragments that cross the boundary. Probes placed very close to a boundary may yield more boundary crossing fragments. Whatever the probe distribution, the count information must do its best to impute. There is no better source of information. 

## 3) Dispersion as a function of the mean -- 

we might be better off completely discarding our overdispersion estimation and shrinkage for now, because it has catastrophic behavior for hybrid capture. However, we could rebuild it as a function of the mean. Suppose we take the counts of all observation nodes (regions and boundaries) and sort them numerically. We could either bin the counts (requires a 'magic number' of bins), use a set of 'nearest neighbor' counts (requires a magic number of nearest neighbors), or even sweep a kernel across the sorted counts (requires selecting a kernel bandwidth, but a gaussian kernel could weight nearer neighbors more than distant neighbors, a nice property). Any of these could be use to create a function returning a 'locally fitted' dispersion as a function of the count. This could be a way to elegantly handle dispersion in the setting of hybrid capture.  

## 4) Imputation variance should be related to the strength of the 'anchoring' observable boundaries that are used for imputation. 

## 5) Count imputation using the gDNA fraction paradigm --- 

this is a promising direction that has not been completely explored. A subset of intron-exon boundaries have both spliced and unspliced fragments crossing the boundary into the adjacent region. For these boundaries, imputation deconvolution can take the form of a fraction (unspliced density / (unspliced density + spliced density)) that is project onto the raw count of the non-observable exon. We are not using this mechanism currently. We are only using the gdna density computed from the count. 

## 6) contained vs crossing count overdisperion modeling -- 
these were previously separated. It is unclear whether there is a role for this anyore. 

--------

In summary, I have several intricate thoughts outlined above regarding the direction we should take. Please assimilate these thoughts. Rewrite a v2 count channel design plan. This time, try to formulate the implementation plan itself.

============================================

Overall thoughts about calibration

Algorithm:
1) Accumulator
2) global RNA FL (spliced)
3) global gDNA FL (seed regions)

4) Strand model

- Fit beta-binomial model (overdispersion)
- Apply strand model to EVERY node

Given sense counts and antisense counts, deconvolve to gDNA counts and RNA counts.

For unstranded data, there should be no 


# Decouple joint model

I am convinced now that we need a robust strand deconvolution module that activates prior to the count deconvolution module. Let's think about it. ## Unstranded -- When we have unstranded data, strand deconvolution does absolutely nothing. It is a no-op. Joint module useless. ## Strong strand-specificity -- Here, we would rather trust our strand deconvolution. The count model tends to make the results worse. We only *need* the count module for ambiguous regions. This is what led to the overdispersion modeling worse to weaken the count model's effect in the joint deconvolution. The count x strand joint was not being weighted appropriately, so a mediocre count model was much strong than an excellent strand model. Titrating this is turning out to be quite challenging! ## Weak strand specificity -- here is probably the only situation where we can benefit from synergy. In real RNA-seq data, "weak" strand specificity is quite rare. Strand-specific data is often 99% strand-specific. Occasionally, I've seen >90%. I'm not enthusiastic about supporting an over-engineered architecture that uses strand to "clean" and uses strand again. It's complicated. It's not elegant. It's almost never needed. It's hard to parameterize. If we decouple strand and count, then we can focus on making strand model work extremely well, and we can focus on the count model working extremely well. The two should rarely conflict (hopefully). 

Open questions for you -- 

1) if we decouple strand and count, how do we parameterize the strand deconvolution model? We want data to "pass through" the strand model as a "no-op" when we have unstranded data. Essentially unmodified (no deconvolution).

Do we use aspects of the [text](strand_clean_robust_deferred.md) document?

We already have the beta-binomial fit. That might be the most robust strand deconvolution. Is that what we should use? 

2) when we have an excellent strand model, how do we make sure the count model doesn't 'ruin' already excellent strand deconvolution predictions? conversely, if we have bad strand predictions.




# Ongoing development of blended strand and count channels

Let's think about the behavior of the new "blending" weight. For each node we have raw counts. After running strand deconvolution we will have deconvolved gDNA + RNA counts. The deconvoluted gDNA count will be less than or equal to the raw count (0 <= deconvolved gDNA count <= total raw counts). Do we "blend" this with the raw counts next? So we have a blending of raw counts and gDNA deconvoluted counts based on the blending weight? Raw counts is always greater than gDNA counts, so the blend is always going to be greater than the gDNA count value. Is that how the behavior works? Does this happen before imputation of non-observable nodes? If this is the algorithm, we can infer a few properties. First, blending raw count with deconvolved count will always be an upward bias to gDNA counts. This may be very minor, especially when strand specificity is high. Second, blending strand-deconvolved counts with raw counts using the blending weight effectively "cleans" the count data. After this blending, we have our "cleaned" count for imputation. This is Phase 3. So is Phase 3 even necessary? It seems like we get Phase 3 "for free" just by running strand deconvolution and blending it with raw counts.

Audit the calibration code. Tell me if my understanding is correct, and if not, what is my misunderstanding? 


# Claude Response

1. The shrinkage target is the count estimate, not the raw count
This is the key correction. The blend is one operation over two gDNA-fraction estimates of the same quantity:


gDNA_blend = w·(g_strand·M) + (1−w)·(g_count·M)     # both ≤ M; M = raw count
g_count is the count module's gDNA estimate (density·eff/M), not the raw total M. So in your toy (strand says 20 gDNA; count imputes ~10): the blend is 0.96·20 + 0.04·10 = 19.6 — it shrinks toward the count's 10 (downward), not toward the raw 100 (up to 22).

Why this matters: "shrink toward the raw count" = shrink toward "all unspliced is gDNA" (fraction 1.0). At low SS (w→0) that collapses to everything-is-gDNA — which is precisely the old cliff that manufactured the gdna_none false positive. The robust fallback at low SS is the count's gDNA estimate (≈0 for a pure-RNA library → no FP), not the raw total. So your shrinkage concept is right (precision-weighted toward a more-robust estimate), but the target must be the count's gDNA estimate. There's also just one blend, not a strand→raw blend followed by a count→strand blend.

2. Strand already deconvolves exons — Phase 1 gives you that; Phase 3 is smaller than it looks
Your core motivation — "use the accurate strand model to deconvolve exons, instead of the count model imputing exons from boundaries" — already happens in Phase 1. Exons of a stranded gene are strand-observable, so the blend computes g_strand for them directly (w=0.96). We are not imputing stranded exons from boundaries; strand deconvolves them in place. So that win is banked.

What Phase 3 (clean→count) actually changes is the count module's minority role:

the (1−w)=4% count component at stranded exons, and
the AMBIG / strand-blind nodes that are pure count.
Its real value is preventing nascent→gDNA mis-calls: under nascent, the raw boundary crossings include stranded nascent, so the count imputation counts nascent as gDNA. Cleaning the anchors fixes that — but bounded by the count's weight. So in stranded data Phase 3 is a modest refinement, not a big accuracy jump. And critically: the dominant residual (unstranded+capture, ~18%) is untouched by Phase 3 — there's no strand to clean with, so the count imputes exons from raw boundaries regardless. That residual is a count-model accuracy problem (the point-5 unspliced-fraction projection / anchor-strength), where strand cannot help. I want to set that expectation honestly.

3. Boundary transport: what it is, and why it can't move before the blend
_transport_boundary_flux runs in assemble_priors, after calibration, on the deconvolved per-region gDNA masses (contained + left + right). For each internal boundary it pools the two sides' gDNA and re-attributes it ∝ density_ratio · boundary_capacity — i.e. the higher-local-gDNA-density side (the probed exon under capture) claims the crossing mass off the depleted intron side. It iterates (the densities update as mass moves), conserving total mass.

So it consumes the deconv output → it cannot run "before blending" (the blend produces its input). It's a post-deconv mass-attribution step, distinct from the count density-imputation step. With Phase 3 the masses entering it are already cleaner, so they compose. But your instinct points somewhere real: the transport's density_ratio·capacity attribution is exactly the smarter logic the exon imputation lacks — the count module currently imputes an exon by averaging its boundary crossing densities, when it should attribute crossing mass to the enriched side the way the transport does. So the right move isn't to relocate the transport; it's to give the count imputation the transport's density-ratio attribution (a point-5 / Phase-4 improvement). That would directly attack the unstranded+capture under-estimate.

4. The weight: the steep curve should emerge, not be a contrived sigmoid
(2κ−1)² is the strand Fisher information per fragment — depth-blind, which is why it caps at 0.96 at ss0.99. The rational way to "trust strand even more at high SS" is to restore what (2κ−1)² dropped: depth and the count's reliability, via inverse-variance/MSE weighting:


w = I_strand / (I_strand + I_count)
   I_strand  ∝ N·(2κ−1)² / (overdispersion)     # ← (2κ−1)² gives rapid low-SS degradation
   I_count   = 1 / (σ_count² + bias_count²)      # ← SMALL where count is unreliable (capture)
High SS, well-covered, capture-biased count → I_strand large, I_count small → w→1 (exactly the flagship behavior you want).
Unstranded → I_strand ∝ (2κ−1)² → 0 → w→0 (rapid degradation — same shape you like).
So the steep-at-high-SS / rapid-degradation curve falls out of the statistics once the count carries an honest variance (Phase 4) — and the key is that under capture the count's variance is large (uncertain imputation), which pushes w→1 without any contrived shape. A sigmoid would be a heuristic stand-in for this; I'd avoid it. The honest count posterior (Phase 4) is the principled "sigmoid." (And note: this reframes Phase 4 as the fix for the +1pp regression and the weight-steepness ask simultaneously.)



# newest prompt

We are going to need to accept this effective length fix. It's a key bug fix. Your metrics and characterization of the benchmarks need to change. Previously you looked at gDNA vs RNA. Now we need to look more carefully at gDNA vs nascent RNA (nRNA) vs mature RNA (mRNA). I suspect that the "leak" looks much worse at the DNA vs RNA level, but that the "leak" is not necessarily worse at the gDNA vs mature RNA level. I suspect that both gDNA and mature RNA could be leaking fragments into the resurrected nascent RNA component. This is not catastrophic. We need to look at our mature RNA false positive rate. If it remains controlled, then we are uncovering the TRUE identifiability challenges that we managed to bury because of this bug. The fact that we uncovered a true identifiability issue is completely different from bugs in the code and shouldn't force us to hold pushing new code that fixes a bug. So tasks -- 1) update your pool level metrics to report gDNA <-> nRNA <-> mRNA as separate pools. Characterize "leak" across the three pools. One pool can leak to either of the other two. 2) Yes, dig into the EM as you proposed to ensure that we don't have other latent bugs related to fragment length discrimination. Our EM has an excellent FL likelihood model that should be helping us.