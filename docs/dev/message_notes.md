Responses:
# 1) Is r ever load-bearing, or is it always a proxy for "a level I could not compare"

With hybrid capture ON, we may improve performance by imputing composition instead of abundance (density e.g. LEVEL). It depends on the object that we are solving and the circumstances.

# 2) Does anything need to flow OUT of an exon?

For simple single-stranded exons, the answer is no, nothing *needs* to flow out of an exon.

However, there are extremely complex loci with overlapping transcripts on opposite strands and numerous isoforms of all kinds. It becomes extremely challenging. I'll give you an example:

Locus Example:

TA+ exons (1000, 2000), (3000, 5000), (10000, 12000)
TB+ exons (3100, 5000), (10000, 12000)
TC+ exons (4000, 5000), (11000, 14000)
TD- exons (7000, 8000)
TE- exons (2500, 20000)
TF+ exons (4000, 5000), (10000, 11000), (15000, 16000), (22000, 25000)

A locus like this stresses the system. Overlapping transcripts in opposite directions, alternate TSS and TES sites. Very complex. These loci are present in the human genome.

# 3) The anchor under nascent RNA

The density model can deconvolve introns. This is effective and a powerful mechanism to solve for nascent RNA when we have unstranded data

# 4) Capture. 

I mostly agree. The physical claim is, "if an exon is either probed or not, nucleic acid sequences overlapping the probe (both gDNA and RNA) are enriched. Enrichment depends on the overlap of the nucleic acid sequence with the probe". There are many probe panels, and most are not published. Most are not easily obtainable. Therefore, Rigel aims to learn the enrichment pattern from the data.

# 5) The 2.5 % boundary-crossing deficit in the queue -- this is likely a bug in the accumulator. It is probably independent of the other problems. It is worth addressing in a dedicated session.

-----------------------------------------

Here is a synopsis of how we can solve single-stranded objects:

# Single-stranded objects

## Intergenic REGIONS: pure gDNA, fixed, static. An individual intergenic REGION measures the local gDNA level.

## Intergenic|exon BOUNDARIES: transcript edges that touch intergenic space. These are TSS or TES. Any boundary crossing fragments are pure gDNA, fixed, static. These measure *local* gDNA levels.

## Intron REGIONS: gDNA + synthetic nascent RNA (unannotated). Two strategies to deconvolute: 1) strand model (best), and 2) density model where we compare the density to the intergenic background density. we have an effective negative binomial model for this

## Intron|exon boundaries

This refers to a BOUNDARY with an INTRON (no exon bits) on one side, and an EXON (exon bits set) on the other side.

The intron|exon boundary is special. It has TWO SIDES: one side facing the intron, and one side facing the exon.

This kind of BOUNDARY has important flags/bits:
- Is this a TSS or TES?
- Do we have a change in RNA strand?
- Is there a splice junction? Which strand is the splice junction on?

We need a lot of careful logic to handle intron|exon boundaries correctly. These boundaries can be partially enriched by hybrid capture if a probe is nearby. So imputing the level from the adjacent intron can be incorrect.

To solve the intron|exon boundary, we can break down the fragments into spliced fragments and unspliced fragments.

To solve for the unspliced fragments, we can impute the COMPOSITION from the intron. If the intron is depleted, then the composition can be poorly define (discrete counts, sparse, can't precisely tell us the composition). It may be challenging to impute nascent RNA and gDNA from a depleted intron.

When solving for the boundary counts, remember, the spliced fragments are excluded. Spliced fragments are pure RNA, they are measured, and they are not needed to solve.

Typically, we spend the vast majority of effort handling the logic at intron|exon boundaries.


## EXON regions (single-stranded)

Solving exon regions depends on what is adjacent. The region will have two boundaries as neighbors. It depends on what the boundaries are.

### intron|exon boundary: The boundary is composed of spliced fragments and unspliced fragments. The gDNA composition is a share of the unspliced fragments. The RNA composition combines the unspliced RNA crossing fragments and the spliced fragments. The Boundary can form a composition and it can be used to impute the exon composition.

### intergenic|exon boundary: and intergenic|exon boundary has only gDNA crossing fragment. Only the gDNA level can be transferred

### exon|exon boundary: typically, we would transfer by LEVEL (not composition) across an exon|exon boundary. These boundaries are typically TSS/TES and occur when strands change.

=========

The synopsis above is just the beginning. It covers the general message propagation for single-stranded objects. Both-stranded objects add complexity but the rules don't necessarily changes. The individual strands can often be handled independently. Deconvolution, however, is 2 degrees of freedom and benefits from having the gDNA landscape prior built from pass-0.

==========

Responses to your questions:

# 1) Capture (the 11.7× row)

- what is the 'level ledger'? I use the term *density* for (counts/bp) units, but the best term is probably 'abundance' for (counts/bp) units. you're using 'level'. We really pushed to constrain our vocabulary and be more precise with words. let's be consistent.

- if we can successfully impute *some* (not all) exons, we can fit a gdna landscape (it is not NPMLE, it is a poisson kernel). A gdna landscape fitted to capture-ON data is bimodal. A correctly fit bimodal gdna landscape is powerful. In the second pass, it pulls gdna estimations towards the correct mode (enriched vs depleted). the lesson is -- we do not have to accurately solve ALL exons during pass-0. We only need enough to train the gDNA landscape. Naturally, the simplest exons to solve are single-stranded exons that are bordered by splice junctions. The splice junctions give us MEASURED pure RNA, and the unspliced fragments crossing the same boundary give us a fairly accurate measurement of gDNA. Together these exon|intron boundaries give us a composition that we can impute onto the exon. This is likely our most promising way to solve a subset of exons; ideally, enough to fit the gdna landscape prior (poisson kernel).

- when message propagation reaches an enriched exon, we have to determine how abundance is interpreted. this should be the job of the message recipient. The message recipient (the destination) knows its population, and knows how it must interpret messages.

- I think there are TWO strategies for message recipients to receive messages. Let's call these the 'ABUNDANCE' strategy and the 'COMPOSITION' strategy.

## ABUNDANCE strategy

- Destination receives message with { gDNA, RNA+, RNA- } abundances (also called densities and levels)
- Destination compares its TOTAL ABUNDANCE against the source node's total abundance (the same thing we do for composition transfer)
- The 'disagreement' can be measured as the ratio of dst_abundance / src_abundance, or log(dst_total_abundance) - log(src_total_abundance). The greater the disagreement, the less trustworthy the message is. 
- If there is a huge disagreement, the destination needs to say, "hmm.. I'm getting a message with a dramatically different set of abundances. I can't trust it"
- This is where the TRANSFER VARIANCE comes into play. The transfer variance is based on the disagreement between source and destination abundances. When there is a huge disagreement, the solution is to dampen the message precision.
- Messages are integrated using a precision-weighted combination, so a high disagreement -> high variance -> low precision -> message dampened -> destination node ignores message.
- This has been the strategy that has worked the best for us in the past.
- We have had a few ways of calculating the transfer variance based on disagreement. The current method is based on building a total density landscape at initialization time. The total density landscape gives us a sense of the global density spectrum and helps us to compute an honest transfer variance between source and destination. I believe this is called the 'NPMLE' at present. It is quite possible that we can use a different strategy. In fact our poisson landscape approach may be well suited for this task as well (we build the poisson landscape from gDNA AFTER the pass-0 solve, but we could build it from total abundance BEFORE the pass-0 solve)

## COMPOSITION strategy

- I outlined how the composition strategy works already, and I believe you've integrated parts of it.
- To use the composition transfer strategy, the **destination transcript population must be the same as the source transcript population**. If we can guarantee that the transcript populations are the same between source and destination, then we reason that the disagreement between total abundance must be due to enrichment/depletion by capture probes and not because of a change in transcript population.
- When we use the composition strategy, we RESCALE the source or destination by the ratio of their total abundances. After rescaling, the abundances become directly comparable.
- We need to carefully think about the transfer variance of the composition model. The precision of the source message must be honest. Scaling an imprecise set of abundances doesn't make them more precise.. the precisions need to be handled correctly. This may require some derivation work.
- The "disagreement variance" concept from the ABUNDANCE strategy does not apply to the composition model, because the message propagates after rescaling such that the total abundances are equal. There is a different derivation and model for the transfer variance under the COMPOSITION strategy. It may be based on the rescaling itself. For example, if you have gDNA = (2 fragments / some length) and a scaling factor of 122X when jumping from a depleted to an enriched node, we get 2x122=(244 fragments / some length). But the precision is still based on 2 fragments. If you think about it, the problem comes from multiplying a discretely sampling integer count by a factor. What we are really doing is upsampling or downsampling. If we have abundances { gDNA = 2 frags, RNA+ = 8 frags, RNA- = 0}, and we have a composition transfer to a destination with { gDNA = 200 frags, RNA+ = 1000 frags, RNA- = 0}. total abundance ratio is 1000 / 10 = 100X. But effectively we are doing binomial upsampling to equalize the abundances, and this produces a very imprecise estimate. I think this is the kind of variance model we need to derive.

## 2) RNA arms at a terminus

Yes, a terminus (TSS or TES) pertains to each strand. We could have a boundary with a terminus on both strands, or a boundary with a terminus on just one strand. 

My suggestion is to let the message recipient determine what to do with the message. The sender can blindly send messages. This is mostly semantics but it helps me think about the problem.

Termini occur at boundaries, and so this situation happens when the message recipient (destination) is a boundary with a terminus. For me it gets quite abstract to discuss it without an example, so here is an example.

Transcript TA+ exons (1000, 10000)
Transcript TB+ exons (5000, 10000)
Transcript TC- exons (500, 20000)

Let's examine the Boundary 'B' at pos=5000. This boundary has both RNA strands active and it is the TSS for TB+. Let's called region `R1` at (1000,5000)

- R1 source sends message to B { gDNA = 1, RNA+ = 10, RNA- = 100 }
- Boundary 'B' has a TSS on the (+) strand. TB+ starts at boundary `B` so none of the crossing fragments can be attributed to TB+.
- `B` receives an incoming message with { gDNA = 1, RNA+ = 10, RNA- = 100 }
- The transcript population at (1000, 5000) is { TA+, TC- }
- The transcript population at boundary `B` is { TA+, TC- }... remember TB+ STARTS at boundary B. It is not part of boundary B.
- **Therefore, transfer can be by composition for R1 -> B**

Now, lets call region R2 (5000, 10000) and continue with forward message propagation.

- B sends message to R2.
- R2 receives message. 
- Transcript population at B is { TA+, TC- }
- Transcript population at R2 is { TA+, TB+, TC- }
- Populations are not the same
- Transfer must be by ABUNDANCE strategy (e.g. "level")

Suppose the abundances at B = { gDNA = 1, RNA+ = 10, RNA- = 100 }, And the total abundance at R2 = 10000. There are two explanations for the discrepancy, one that TB+ has very high abundance, or two that there is a capture probe on TB+ that enriches all of the components. At pass-0 we don't know.

So for the B -> R2 message propagation, R2 receives the message. What should R2 do?
- R2 thinks, "wow, there is a huge discrepancy here. i have TB+ and this message doesn't have TB+. This could either be TB+ or it could be hybrid capture, or BOTH"
- R2 has total 10000, whereas B2 message { gDNA = 1, RNA+ = 10, RNA- = 100 } total abundance 111.
- What should we do here?
- By my abundance strategy, we would transfer the abundances, but this would leave the destination R2 with (10000-111) = 9889 abundance unaccounted for. There is not enough information to assign this abundance to any transcript or gDNA.

Ultimately, this example will prove unsolvable! The reason is because all of the transcripts are single-exon transcripts. We have no measured RNA anywhere in the system. We need splice junctions to give us spliced RNA to solve this example, because it will give us measurements for the transcripts that otherwise have no measurements.

Even with a correct gdna landscape prior, we still won't be able to solve, because we need a way to estimate the gDNA level in order to project onto the prior landscape.

- If we have STRAND SPECIFIC data, we CAN solve this system with a gDNA prior. 

I hope this at least provides insight into the mechanisms involved.

-----

## 3) The rescale's totals

Again I suggest that the message SOURCE just sends the message without any regard to what happens to it. The destination needs to decide how to interpret the message. My example was SPLICE-OUT, correct? 

So we have this example:

EXON REGION R1 <-> EXON|INTRON BOUNDARY B <-> INTRON REGION R2

- R1 sends message to B (R1 -> B).
- B receives message, sees that we can use composition transfer
- B rescales by TOTAL ABUNDANCE ratio (R1 / B), still no belief, just measured total abundance **** SEE ASIDE ON THIS ****
- AFter the rescale, R1 and B are on equal total abundance scale. Transfer message R1 -> B (use the upsampling vs downsampling concept, transfer variance/precision)

- Now, B has ONE MESSAGE from R1. It still has no belief.
- During message propagation, B has no belief, so it RELAYS the message, unmodified (the intron HAS an initial belief, so the message will likely dampen or die there)
- B still can't solve because it needs its other message (from the intron).
- The intron HAS a belief (it deconvolves into gDNA + synthetic nascent RNA).

FINALLY, at solve time, B has two messages, ONE from R1, and ONE from R2.
- R1's message has been rescaled
- R2's message must also be rescaled (also composition transfer compatible)
- B has no belief of its own
- So the average of the precision-weighted messages becomes the B belief.

## 4) Rung 3

- A single-stranded INTRON can deconvolve by using the density model (this rests of the assumption that intergenic regions measure background gDNA, and introns measure background gDNA plus synthetic nascent RNA).
- Introns can solve this way at initialization time, so that they have a belief to propagate at pass-0
- Note that with strand-specific data, ANY single-stranded nodes (including introns) can solve and express a belief as well

## 5) Naming - 'CurrencyPolicy' makes no sense to me. Right now the name doesn't matter. We will eventually replace 'RelayPolicy' once we finish development.

=============

Yes, as we think of additional transcript cases, we should add them to the test chromosome. This will evolve and become a gold standard set of 'goldens' that can be added as a type of regression test for the tool. We must keep adding test cases whenever we think of a new situation that could stress the tool or expose a weakness.

================

I have given you a lot of information in this turn. I want you to be meticulous in your considerations and in your design. Organize your thoughts carefully. Plan your implementation carefully. Test one new concept at a time. Be tenacious and methodical

============================

# Total density (total abundance)

Knowing the total abundance (counts/bp) of each node is critical, because we compute enrichment ratios based on total abundance.

If we store 'mass' (fractional counts), we must convert it to an abundance by dividing by the effective length. However, the effective length depends on the composition. gDNA and RNA have different FL distributions.

If gDNA has mean frag len = 100 and RNA has mean frag len = 200, we have a big discrepancy. If we have 100 total counts in a 500bp region, then ROUGHLY:
- If 100% gDNA, then 100 counts / (500-100) = 0.25
- If 100% RNA, then 100 counts / (500-200) = 0.33

This problem is greatest at boundaries, because boundary-crossing fragments directly depend on fragment length.


## Goal: compute total abundance at accumulator time

During the accumulator step, we deposit fragments with known length into regions and boundaries. We deposit a few values:
- counts (incidence), used for statistical power
- mass (conserved), when a fragment overlaps many boundaries, it deposits mass that sums to exactly 1.0
- abundance!?

We have gone through several iterations of depositing abundance. For the message propagation to work, we need to compute enrichment ratios using total abundance (gDNA and RNA combined).

For (*most*) boundaries, given a fragment length 'w', the abundance (counts/bp) is `(1.0 / (w-1))`. There are some caveats, such as when a boundary is very near the end of a transcript, that constrain the start positions of a fragment with length 'w'. However, we already went through **MANY** iterations of trying to model this and proved that it is a circular problem -- depositing correct abundance for the corner cases (where the fragment length is constrained by transcript ends) REQUIRES knowing that the fragment is RNA (gDNA has no such constraints) and knowing which transcript isoform the RNA belongs to. That is the whole point of the tool. So depositing boundary abundance `1.0/(1-w)` is roughly correct with some known limitations.

For regions with length `L`, the effective length for fragment deposition is `L - w`. There is a derivation (we derived this) for total abundance for depositing contained fragments in REGIONS. It may be useful to find this in the docs and perhaps the delete aspects of the git repo to find it, or re-derive it.

My question to you is -- how are we trusting our enrichment ratios now? I believe our simulated scenarios have equal gDNA FL and RNA FL so this problem doesn't expose itself. IT's a major problem when FL is highly different!

This might be something we have to address now.

===========


I like the idea of delivering both claims with honest variances. This is an elegant solution. It lets the data decide how to interpret messages. 

Interestingly, I believe the absolute abundance transfer and the composition transfer are on two ends of a spectrum. The two strategies are deeply connected somehow, and my intuition tells me that there is actually a knob that connects them. Absolute abundance transfer is possible when *total* abundance is uniform and involves addition/subtraction. Composition transfer is required when *total* abundance is non-uniform (enrichment/depletion) and involves multiplication/division. I think, it may never be exactly one or the other, sometimes it is both. And you seem to be honing in on that idea as well. At least, your first idea is a 'switch' that chooses between two message propagation strategies. I wonder if the evolution of the idea will take us somewhere even more elegant.

========================

You question about 'total-abundance landscape built before pass-0'. Right now, I don't think we need the total abundance landscape because the local disagreement `σ² = (log(T_dst/T_src))²` is a reasonable place to start. The problem with local disagreement is that it has no grounding. There is no concept of "who is right? who is wrong?" Projecting onto a total density landscape projects an observed abundance at a single node onto the global landscape of abundances. The assumption is that there modes (enriched mode, depleted mode, etc) under hybrid capture and observed abundances are sampled from those modes. Projecting onto the global abundance landscape finds the nearest mode, so instead of computing variance from the local abundance disagreement, we compute variance based on the projections of the local abundances onto the global landscape. For example, if we have sparse counts and two observations are just sampled differently from the same mean, they may have high disagreement. Projecting onto global landscape helps us robustify the individual local observations.

I'll reiterate, I don't think it's essential for the first iteration of the policy. You can implement and we can refine later.

=============

A few things:
1) I'd like to compare the NPMLE against the new abundance landscape. I want to see what the landscape plots look like head to head. Which is more accurate?
2) I look at the abundance landsape plots. They are kind of "spiky" with more modes that I would expect. This is not necessarily wrong, but setting the bandwidth of the poisson kernel seems like a parameter worth investigating further. How do we know our bandwidth is optimal?
3) What can we use our abundance landscape for in the pass-0 solve? Remember the role of the NPMLE --- the NPMLE was used in pass-0 to estimate the "transfer variance" for message propagation. Please audit this behavior. How does the transfer variance estimation work? How does it work with NPMLE versus with our new abundance landscape? Will this be better or more accurate?



====================================

# New calibration plan

# Calibration First Pass

## 0. Categorize regions and boundaries

The first pass needs a SUBSTRATE to train a gDNA prior. The substrate must include exons. We need to identify regions and boundaries (including exons) that can be solved and used for training.

We can identify the substrate purely based on structural properties of the transcriptome (at index creation time).

The substrate includes:
- Intergenic regions and intergenic-exon boundaries
- Single-stranded introns
- Single-stranded intron-exon boundaries

The EXON substrate must be carefully defined, and we must ensure that we correctly define our 'solvable' substrate for the pass-0 phase.

EXONS that are universally solvable include:
- single-stranded
- at least one side borders an exon|intron boundary with a splice junction. This boundary must have no contiguous exon crossing (no exon|exon boundary component).
- ok if the exon has two exon|intron splice junction boundaries

### Example 1 of solvable exons:

Here are two transcripts:
- Transcript TA+ exons (1000, 2000), (10000, 11000)
- Transcript TB+ exons (1000, 3000), (10000, 11000)

Here, the EXON (10000, 11000) qualifies as solvable, because it is single-stranded and has the splice junction exon|intron boundary at position 10000.

### Example 2:

Given transcripts:
- Transcript TA+ exons (1000, 2000),(4000,5000),(10000,11000)
- Transcript TB- exons (500, 3000), (10500,15000)

There are no solvable exons here. Everything is both-stranded (AMBIG)


### Example 3:

Given transcripts:
- Transcript TA+ exons (1000, 2000),(4000,5000),(10000,11000)
- Transcript TB+ exons (500, 3000), (10500,15000)

Solvable exon regions include:
- (2000, 3000) region bordered by exon|intron sj boundary at 3000
- (4000, 5000) single-stranded, TWO sj exon|intron boundaries on both sides 4000, and 5000.


## 1. Intergenic initialization

The first stage of calibration is initializing the structurally fixed regions:

- Intergenic regions: 100% gDNA, variance = 0, fixed
- Intergenic|exon boundaries: 100% gDNA, variance =0, fixed

We then *model* the intergenic gDNA distribution from intergenic REGIONS (not boundaries). Then move to intron deconvolution.


## 2. Intron deconvolution

The second stage of calibration is intron deconvolution. 

*Eligible introns must be single-stranded*

Intron deconvolution requires:
1) strand model (uninformative when ss=0.5)
2) density model (intergenic region distribution)

We apply our abundance model + our strand model to deconvolute introns. The solve is message free (no outside information needed). 

**What should the prior be?** The question is, "do we need a prior here?"

- predicted fraction of gdna = (intergenic gdna density) / (total density of intron). this can be the prior, right? if the intergenic gdna density is greater than the intron total density, then we clamp to 100% gdna as the composition prior.

**The intron solver needs to integrate the strand model and the density model** The abundance model uses negative binomial model to provide a solution based on total abundance (fraction gdna, fraction rna) and the strand model can independently provide a solution. The two results must be integrated using an honest precision-weighted combine.


## 3. Intron|exon boundary deconvolution

**The goal of this step is to solve intron|exon boundaries that are adjacent to solved introns**

- Every solvable intron now has a belief established from the strand model and the density model
- The next step is to solve the intron|exon boundaries that are adjacent to solved introns.
- **This step is technically a message passing step**

We "fan out" from the solved introns in step 2 (above) to their neighbors.

The layout looks like this:
<exon> | **<intron|exon boundary>** | <solved-intron> | **<intron|exon boundary>** | <exon>

**How do we solve for the intron|exon boundaries?**

### Composition transfer 

We borrow the concept of 'composition' transfer. If the intron belief is (85% gdna, 15% rna), we reason that the adjacent intron|exon boundary must share the same composition belief, because the fragments crossing from the intron to the intron|exon boundary are distinctly intronic and share the composition.

Composition transfer still needs to be treated statistically and assigned honest precision. Why? Because ultimately, we have multiple channels of information! The intron|exon boundaries could independently solve with the strand model, without requiring any information from the neighboring introns.

### Example: count discrepancy

- suppose the intron has 3 fragments. deconvolution with the strand model and density model must ultimately assign fragments (a fragment is a discrete observation, and must be either gDNA or RNA). So we only have a few possibilities for the 3 fragments (3/3 = 100% gdna), (2/3 = 66% gDNA), (1/3 = 33% gDNA), or (0/3 = 0% gdna). One degree of freedom.

- suppose the deconvolution predicts (33% gDNA, 66% RNA) because the gDNA background density low.

Now, we look at the neighbors. Suppose the adjacent intron|exon boundary has 300 unspliced crossing fragments (100X). How do we impute?

We are imputing a (33% gdna, 66% rna) composition derived from only 3 fragments to a 300 fragment exon|intron boundary. The arithmetic would set this to 100 fragments gdna, 200 RNA. But what about the precision?

#### Composition transfer proposed derivation

To carry honest precision during upsampling from 3 fragments to 300, you must model the composition as a **posterior probability distribution** over proportions rather than a fixed point estimate, using a **Beta-Binomial framework** or **Law of Total Variance**.

**1. Beta-Binomial Predictive Model**
Instead of assuming the gDNA fraction $p$ is fixed at $1/3$, treat $p$ as a random variable derived from the 3 observed fragments using a prior distribution $\text{Beta}(\alpha_0, \beta_0)$ (e.g., uninformative flat prior $\alpha_0=1, \beta_0=1$):

* **Posterior for $p$:** Given $k_1 = 1$ gDNA out of $n_1 = 3$ fragments:

$$p \mid k_1 \sim \text{Beta}(1 + \alpha_0,\, 2 + \beta_0) = \text{Beta}(2, 3)$$

* **Predictive Distribution at Boundary ($N_2 = 300$):** Impute the gDNA count $K_2$ using the posterior predictive distribution:

$$K_2 \sim \text{Beta-Binomial}(N_2 = 300,\, \alpha = 2,\, \beta = 3)$$

**2. Quantifying Uncertainty (Variance Propagation)**
By the **Law of Total Variance**, the total variance of your imputed boundary count $K_2$ splits into two components:

$$\text{Var}(K_2) = \underbrace{N_2 \cdot \mathbb{E}[p(1-p)]}_{\text{Binomial Sampling Noise}} + \underbrace{N_2^2 \cdot \text{Var}(p)}_{\text{Epistemic Uncertainty from } n=3}$$

* **Naive Point Estimate ($p = 0.4$ fixed):** $\text{Var}(K_2) = 300 \times 0.4 \times 0.6 = 72 \quad (\sigma \approx 8.5 \text{ fragments})$
* **Beta-Binomial Predictive Model:** $\text{Var}(K_2) = 3,660 \quad (\sigma \approx 60.5 \text{ fragments})$

Because the epistemic uncertainty scales quadratically ($N_2^2$), the low count at $n=3$ dominates the calculation. Instead of precise $(120 \text{ gDNA}, 180 \text{ RNA})$ counts, the honest 95% Credible Interval for gDNA spans roughly **18 to 228 fragments**.

**Practical Implementation Strategies**

* **Precision Weighting:** If fitting a downstream model across junctions, weight the boundary observation by the inverse variance $w_i = 1 / \text{Var}(K_2)$ so the system automatically discounts this highly imprecise upsampled boundary.


### composition transfer solve overview

We solve for an intron|exon boundary.

Inputs:
- adjacent neighbor intron composition. the raw counts (gdna fragments, rna fragments) give us our statistical power
- intron|exon boundary total counts, counts per strand

Procedure:
- up-sample / down-sample from the intron to the intron|exon boundary to match total counts.
- the precision model described above could be used to do the upsampling/downsampling while maintaining honest precision based on discrete counts (each count comes from either gdna or rna).
- precision-weighted combine between the composition transfer solution (upsampled/downsampled composition with precision)
- solve the intron|exon boundary using its own strand model information (uninformative if we have unstranded data) and the composition transfer from the neighboring intron




## 4. EXON region deconvolution

This is stage 4. To reach this stage, we have already initialized intergenic regions. Then we solved for intronic regions (abundance + strand). Then did a "one-hop" imputation between each solved intron region and the adjacent (intron|exon boundary).

Now, we have enough information to solve a specific subset of EXON regions. The exon subset is described above (part of stage 0 is identifying eligible exons).

Principles:
- Every solvable exon has two boundaries (left, right).
- At least one of the boundaries must be a splice junction AND be an intron-only boundary (no exon-exon fragments crossing contiguously)
- For the eligible boundaries, we can fully interpret the unspliced crossing fragments as gDNA or nascent RNA, AND we have access to SPLICED fragments (pure RNA)


### Solvable exons must be matched with solved intron|exon boundaries

The key step here is that we must identify cases where solvable exons can be matched to solvable introns.

Here is an example that demonstrates the matching:

TA+ exons (1000, 2000), (10000, 11000)
TB+ exons (4000, 5000), (15000, 16000)

In stage 2 we do intron deconvolution:
- intron region (2000, 4000) is solved
- intron region (5000, 10000) is solved
- intron region (11000, 15000) is solved

How do introns match to solvable exons?

- intron region (2000, 4000) matches with exon (1000, 2000) because the (intron|exon) boundary at position 2000 is a splice junction boundary
- intron region (2000, 4000) DOES NOT match with exon (4000, 5000) because the (intron|exon) boundary at position 4000 IS NOT A SPLICE JUNCTION (it is a transcript start site or TSS for transcript TB+)

- intron region (5000, 10000) matches with exon(4000, 5000) AND matches with exon(10000, 11000).

- intron region (11000, 15000) matches ONLY with exon(15000, 16000) because the left side exons (10000,11000) is a transcript end site (NOT A SJ)


### EXON solve procedure: composition transfer again

**this is also a message propagation step**

We will demonstrate the composition transfer between ONE boundary and the EXON region itself. Each solvable EXON has two boundaries. If both the left and right boundaries are eligible, then the composition transfer should be performed for both boundaries. The final solve is the precision-weighted solve from the boundaries AND the strand model (uninformative for unstranded).

Composition transfer procedure:
- exon|intron boundary has already been solved, and has  (gdna counts, unspliced rna counts).
- exon|intron boundary also has SPLICE JUNCTION counts (these are pure RNA), spliced rna counts. The spliced counts are measured with count precision. This is different from the unspliced count precision.

#### There must be a precision-weighted combine within the boundary unspliced + spliced before composition transfer

The boundary has unspliced (gdna counts, rna counts) with precision determined from the previous stage solve, from the intron <-> intron|exon boundary deconvolution and the strand model. 

Before we can perform composition transfer to the neighboring exon, we must do a precision-weighted combine between the unspliced fragments and the spliced fragments. The splice fragments have 100% RNA at directly measured count precision.

So we are merging (unspliced gDNA counts, unspliced RNA counts) at one precision with the (spliced RNA counts) at measured count precision ---> new composition (gDNA, unspliced RNA + spliced RNA) with an integrated precision for transfer.

Conceptually, the spliced fragments have count precision and will drive the precision of the message. But, think about if there are very low spliced fragments and high unspliced. Then the precision is driven by the unspliced fragments. I don't think the precisions are additive. The higher precision value is pulled down by the lower one.


### EXON composition transfer

The incoming data from the boundary has integrated precision from (gdna counts, unspliced rna counts + spliced rna counts). This gives us our new composition. This is what must be transfered.

We then perform "binomial upsampling/downsampling" so that the total counts are equal. The variance of the resampling process becomes the variance of the transfer. This implies that the transfer has a "price" associated with it. Whatever this is becomes the COST of imputation.

So the composition of the exon|intron boundary (gdna fragments, unspliced rna + spliced rna) is transferred. The composition itself has a precision, and then the TRANSFER itself involves an imputation which we model by upsampling/downsampling in count space, which adds a variance or cost of the imputation process itself.

We must do this modeling for each of the eligible boundaries of the EXON region.

### EXON region solve

Finally, we have inputs:
- strand model (internal solve of the exon)
- left boundary (if eligible) composition with precision
- right boundary (if eligible) composition with precision

We perform the precision-weighted integrated solve of all three. *this is a one-hop message propagation event*


===================

After the above, we now have a set of solved regions and boundaries.

Pass-0 thus ends here. Any remaining regions and boundaries are left unsolved. No effort is made during pass-0 to solve them.

Next, we fit the gdna abundance landscape on the solved regions and boundaries ONLY. 

====================

# Open questions

## reference prior?

- when running the reference 

## composition representation: counts or abundance?

- we have two metrics: raw counts (discrete observations) and abundance (counts per base)
- how is composition represented? if composition a ratio of abundances? or is it a partitioning of the raw counts?
- the write-up above assumes composition is a partitioning of the counts.
- if composition refers to the partioning of total abundance into gdna abundance and rna abundance, some of the derivation may (or may not change)
- the relationship between abundance and counts depends on the relative fragment length distributions of gdna and rna


















































