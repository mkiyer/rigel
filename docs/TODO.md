# TODO


## BAM scan pileup

I have an idea that is an offshoot of your "5. Coverage-shape coherence priors" -- during the 


### density model

We are working on calibration and have a roadmap (attached). We need to focus on building a density model that handles both hybrid capture data and non-hybrid capture data. The concept of the problem is well formulated. Goals of calibration are to estimate the amount of gDNA (numerator) and the 'exposure' (denominator). Estimating gDNA requires clever use of 'contained' and 'boundary' fragments in order to respect exposure differences. Hybrid capture changes the exposure itself. Keeping contained and boundary counts separate is a wonderful way to isolate hybrid capture from non-capture behavior. The goal of the density model is to be able to estimate gDNA counts given the data in a region R (a region has spliced pos/neg counts, unspliced boundary pos/neg counts, and contained counts). We need to define our 'training' regions from which we can build models. Our models need to predict gDNA counts given the region size (region end - region start), the boundary flux going into the region, and 'global' gDNA evidence. We can model 'global' gDNA using intergenic and intronic regions, this model must take into account hybrid capture data because intergenic/intronic regions are depleted relative to exons. Hence the need to model gdna counts as well as apparent exposure. 


### gdna density model

Given a region, estimate gDNA counts (with a user specified upper confidence e.g. 95% or 99%)

Model should predict gDNA counts in a region ~ region raw length, region exposure weight, boundary flux gdna evidence, global gdna evidence, and upper tail conf level (e.g. 95% or 99%)


### intergenic and intronic 'contained' evidence

Intergenic regions and intronic regions can be assumed to be gDNA dominant. The 'contained' region counts within intronic and intergenic regions are typically NOT captured by hybrid capture probes. They are typically part of hybrid capture exposure and must be weighted accordingly (we can let the data drive exposure weights).

We need to build a model from all intergenic/intronic regions using the contained fragment evidence


### intergenic and intronic 'boundary' evidence

This is reliable gDNA evidence that is usable in hybrid capture and non hybrid capture cases. With hybrid capture, the boundary flux will be enriched. When we don't have hybrid capture, boundary flux data will be sparse (but still usable).

For each intergenic/intronic region, we have two boundaries (left and right). Each boundary has fractional counts or flux. These can be assumed to be gDNA counts.

Boundary crossing fragments will have the FL distribution of gDNA, so a single boundary has the length of gDNA fragment length distribution. However the boundary exposure also depends on the size of the region, so the gDNA FL distribution must be clipped to region size and integrated/marginalized using our existing algorithmics and functions.

We need to build a model from ALL exon-intron and exon-intergenic boundaries. We model:

boundary crossing flux ~ boundary eff length (gDNA FL distribution clipped by region size)

### effect of hybrid capture

Hybrid capture rna-seq targets transcripts (except for panels that have negative control probes, which are typically a very small number). The challenge here is that hybrid capture uses probes to target certain regions. It is never the case that all exons are targeted (that is exceedingly hard) and so the boundary flux gdna density model will be expected to be a bimodal combination of off-target boundaries with roughly the intergenic/intronic global off-target gdna density and 'on-target' boundaries which will start to approximate the on-target exon density. 

Exon regions can thus be divided into expressed x targeted
1) expressed AND targeted (on-target)
2) expressed AND off-target
3) not expressed, targeted
4) not expressed, not targeted

We can start by building 1) the global intergenic/intronic 'contained' gdna density model and 2) the boundary intron-exon/intergenic-exon boundary crossing model.

We can identify exon group 4 (not expressed, not targeted) by comparing the 'contained' exon counts against the intergenic/intronic gDNA model. "Could we have observed this many contained exon counts by sampling from the intergenic/intronic gdna distribution?". If yes, then we can say the exon is 'not expressed, not targeted'. 

In a 2nd pass these off-target not-expressed exon can be INCLUDED in the global gDNA density model.

Unsupervised learning of hybrid capture means that we still need to partition boundaries into 'on-target' and 'off-target' boundaries.

We should then partition boundaries into 'on-target' and 'off-target'. An easy way is just to compare a boundary density against the global intergenic/intronic gdna density. "Could we have observed this many boundary counts by sampling from the intergenic/intronic gdna distribution"? 

We can then build an on-target boundary flux model from "on-target" boundaries. 

We then have a rough estimate of hybrid capture enrichment by comparing "on-target" boundaries against "off-target" boundaries. 

Finally, we can identify exon group 3 (not expressed, targeted) using our on-target boundary model "if we impute exon counts from boundary counts, could we have observed this many exon counts from the on-target boundary density"?

If we want to, we can then incorporate these "on-target, not expressed" exons into a more robust on-target gdna density model.

We need a simple, elegant procedure to accomplish these tasks. 

Finally, we need exposure weights. In otherwords, we need to factor in the enrichment between 'on-target' and 'off-target' regions into the bayesian prior for EM

OR we need a brilliant idea for handling the denominator for the gDNA component. Our EM normalizes transcripts by effective length and we need some sort of normalization for the gdna component as well, or a competing solution.


Here is my current solution:

This is exactly the kind of problem where principled Bayesian statistics shines. Your instinct to use the boundary flux as independent evidence to predict contained density is the key. It allows us to build a mathematically rigorous model that completely avoids hard classification (e.g., binning exons into groups 1-4) or manual heuristics.

The fundamental invariant we can exploit is **geometry**: for any region $r$, the expected number of unspliced fragments is strictly proportional to its effective length opportunity ($L_{eff}$).

By treating the regional density $\lambda_r$ as a continuous latent variable, we can solve this using an **Empirical Bayes Shrinkage Model (Poisson-Gamma Conjugacy)**. This provides a single, elegant framework that dynamically adapts to hybrid capture targeting.

### The Formulation: Empirical Bayes Shrinkage Density Model

#### 1. The Global Background Prior (Empirical Bayes)
First, we use pure intergenic and deep intronic regions to learn the global "off-target" background distribution.
Let $N_r = N_{contained} + N_{boundary}$ and $L_r = L_{c} + L_{b}$.
We fit a global Gamma prior to the density of these background regions using maximum likelihood:
$$ \lambda \sim \text{Gamma}(\alpha_0, \beta_0) $$
- **$\mathbb{E}[\lambda_{off}] = \alpha_0 / \beta_0$** represents the base gDNA density across the whole genome.
- In hybrid capture, $\alpha_0/\beta_0$ will naturally be very low due to off-target depletion.

#### 2. Local Exposure Inference (The Shrinkage Posterior)
For *any* fine region $r$ (including exons), the boundary flux $N_b$ acts as our independent targeting evidence. 
Since $N_b \sim \text{Poisson}(\lambda_r L_b)$, applying Bayes rule updates our prior, yielding the **exact local density posterior**:
$$ \lambda_r | N_b \sim \text{Gamma}(\alpha_0 + N_b, \beta_0 + L_b) $$

This is where the magic happens:
- **Off-target regions (Sparsity):** $N_b \approx 0$. Because local evidence is weak, density smoothly shrinks to the global background ($\alpha_0 / \beta_0$).
- **On-target regions (Enrichment):** $N_b \gg \alpha_0$. The prior is overwhelmed, and the density confidently snaps to the empirical local exposure ($N_b / L_b$).

#### 3. Predicting Contained gDNA (The Numerator & Upper Tail)
We marginalize this Gamma posterior over the containment opportunity $L_c$ to get the predictive distribution for contained gDNA. Mathematically, Poisson-Gamma conjugacy resolves cleanly to a **Negative Binomial distribution**:
$$ N_{c}^{\text{gDNA}} \sim \text{NegativeBinomial}\left(r = \alpha_0 + N_b,\ p = \frac{\beta_0 + L_b}{\beta_0 + L_b + L_c}\right) $$

- **Expected gDNA:** $\mathbb{E}[N_c^{\text{gDNA}}] = \frac{\alpha_0 + N_b}{\beta_0 + L_b} L_c$
- **The Conservative Upper Bound:** We simply query the exact quantile function of this Negative Binomial distribution (e.g., `scipy.stats.nbinom.ppf(0.99, r, p)`). This elegantly and naturally propagates both *sampling noise* and *exposure uncertainty*!

#### 4. The Exposure Denominator (For EM)
You asked for a brilliant idea for handling the EM normalization denominator for the gDNA component. The EM needs an exposure weight $A_r$. 
This weight dynamically falls out as the ratio of expected local density to global background density:
$$ A_r = \frac{\mathbb{E}[\lambda_r]}{\mathbb{E}[\lambda_{off}]} = \left( \frac{\alpha_0 + N_b}{\beta_0 + L_b} \right) \frac{\beta_0}{\alpha_0} $$
This perfectly represents the true fold-enrichment (or depletion) of every individual region. 

---

### How this resolves your 4 Exon Groups without heuristics
Because this operates on a continuous Gamma posterior, it eliminates the need to do "2nd passes" or hard, brittle partitioning of targets.
1. **Expressed + Targeted:** Huge $N_b \rightarrow$ Posterior $\lambda$ is high $\rightarrow$ Predicts high gDNA. But the true $N_c$ (containing mature RNA) will still vastly exceed the 99th percentile bound.
2. **Expressed + Off-target:** Low $N_b \rightarrow$ Posterior $\lambda$ shrinks to low background. True $N_c$ is high. Massive excess goes to RNA.
3. **Not expressed + Targeted:** Huge $N_b \rightarrow$ Posterior $\lambda$ is high. True $N_c$ perfectly matches the Negative Binomial prediction. ZERO excess for RNA.
4. **Not expressed + Off-target:** Low $N_b \rightarrow$ Shrinks to background. $N_c$ is low. ZERO excess for RNA.

### Robustness to Nascent RNA
Wait, if a region has extreme *nascent RNA* expression, won't nascent RNA cross the boundary and inflate $N_b$? 
**Yes—and that is actually a feature.** Because nascent RNA shares the exact same geometric invariant ($L_c / L_b$) of unspliced gDNA, inflating $\lambda_r$ causes the model to conservatively subtract *both* gDNA and nascent RNA together automatically as the "unspliced background." The remainder ($N_c - k_{\text{max}}$) will securely and purely represent **mature RNA**, fully protecting your EM transcript abundance vectors from intronic/nascent noise.

This framework gives you everything you asked for: mathematically pristine, no magic numbers, robust to hybrid capture, and a native way to surface continuous exposure weights to the EM.



## Phase 4


### rna-seq data types

we essential have four types of rna-seq data: 
1) hybrid capture versus no hybrid capture and 
2) unstranded vs strand-specific -> this gives us four cases that we need to handle: 
a) unstranded, no hybrid capture, 
b) unstranded, hybrid capture (hardest), 
c) strand-specific no hybrid capture (easiest), and 
d) strand-specific hybrid capture (easy). 

Our biggest challenge is designing a set of solutions that we will work each one of the four cases (actually, we can omit unstranded hybrid capture for now because it is rare in real data). So we need to think through the special cases here.

#### no hybrid capture

Your above plan will likely work well for this case. The contained and boundary crossing fragments can be combined into a unified region fragment count and modeled that way

#### hybrid capture

with hybrid capture data, intergenic and intronic regions are depleted. certain exons are highly enriched. other exons are depleted as well (capture only targets certain transcripts, not all transcripts).

for hybrid capture, combining boundary-crossing with contained fragments can dilute our signal substantially, because we combined intergenic (depleted) and intronic (depleted) with exon-intergenic (potentially enriched) and exon-intron (potentially enriched). 

we need an alternative solution for hybrid capture.

I'm happy to start with our most important use case, develop it to fruition, and expand from there. But I don't want our architecture to become too specific so that we can't generalize to the other types of rna-seq data. I want rigel to be able to handle all types of rna-seq

Can you help me with the overarching plan here?


#### strand-specific hybrid capture

this is the most important mode and the most pressing. this is the one we need to start with.

we need to handle all types of strand-specific data because some library prep techniques create different patterns of read 1 and read 2 (sense/antisense vs antisense/sense)

strand-specific data is orthogonal evidence for gDNA vs RNA. we have a global strand model that predicts the strand-specificity of the library. SS = strand specificity which ranges from 0.0 to 1.0. 
- SS = 0.0 means read 1 is 'sense' and read 2 is 'antisense'
- SS = 0.5 is unstranded data (no strand specific information)
- SS = 1.0 means read 1 is 'antisense' and read 2 is 'sense' (most common)

We can easily measure the SS variable from the spliced RNA reads in the library and have an extremely reliable measure of SS. for the sake of describing the approach I am going to refer to SS > 0.99 as highly strand-specific (this is what we commonly see). 

We need to build a gDNA strand balance model for strand-specific data. We need 'training data' to train the model. Training data must come from gDNA.

We can use the following sources of training:
- intergenic (depleted in hybrid capture, sparse)
- intronic (depleted in hybrid capture, sparse)
- exon-intergenic (some enriched in hybrid capture)
- exon-intron (some will be enriched in hybrid capture)
- ***exons with no RNA expression*** how do we assume that an exon has no RNA expression? if we know we have strand-specific data, we can immediately assume that if the number of antisense counts is greater than the numebr of sense counts, the exon is not expressed.. (there is no RNA). Or more generalized, we can come up with a way of use the antisense vs sense breakdown in the exon to determine if the exon is 'expressed' (has RNA) or 'not expressed' (no RNA). we have spliced fragment counts for exons too, so that can be another piece of evidence.

Other key aspects.
- training regions cannot have 'ambiguous' strand. They must have either 'pos' or 'neg' strand.

If we have our set of training regions, then we want to figure out how overdispersed gDNA fragments with respect to strand. gDNA is biological double-stranded and has a concrete fixed mean of 0.50. That is fact. But the overdispersion must be estimated.

We want our model to be able to deconvolute exons that have mixture of gDNA and RNA fragments. We are deconvoluting just the unspliced reads (the spliced reads are pure RNA already).  We want to build a model that can answer the deconvolution question:

"Given that we see `n_ss` counts on the sense strand and `n_as` counts on the antisense strand in this region, how many of the counts can be attributed to gDNA? How many to RNA? Then the next level question is, with 95% (or 99%, or 99.9%) probability, how many counts can be attributed to gDNA?

I am thinking that the variance is dependent on the mean as well, so I am not sure we should model this by strand fraction alone. It's almost like we want to plot log(n_antisense) (x axis) vs log(n_sense) (y axis) and do regression, but that is a simplistic solution. In the past I was doing quantile regression, because I wanted to estimate gDNA with a 95% or 99% confidence. I am not sure that is the best way.


If we iterate through every training regions, we can gather the data for our model. 

- total fragments = contained + boundary left + boundary right
- of the total fragments, how many are 'sense' (same strand as exon or intron) versus 'antisense' (opposite strand as exon or intron)
- total frags = sense frags + antisense frags

can you help me develop this modeling procedure further?

The framing is: Given a region *what is the minimum number of RNA molecules I can assert with 95% confidence?*


Here is some initial work on this. I want you to start from scratch though and derive from first principles.

My intuition is that we can likely model the RNA strand distribution as Binomial safely (we likely do not need beta binomial for RNA). We can test this by iterating through regions and using the spliced fractional counts (pos and neg strand) and estimate the RNA dispersion


### 1. Conservative DNA subtraction (intuitive)

Rather than subtracting the *expected* DNA sense counts, subtract the **upper tail** of what DNA could plausibly produce.

Ask: given `n_DNA` DNA molecules, what is the 99th percentile of sense-strand counts under the DNA beta-binomial?

> k_DNA_max = BetaBin\_quantile(0.99; n\_DNA, μ=0.5, φ\_D)

Then your conservative RNA sense count is:

> k_RNA_conservative = max(0, k_observed − k_DNA_max)

And conservative n_RNA estimate:

> n_RNA_conservative = k_RNA_conservative / 0.99

You're essentially giving DNA every benefit of the doubt — assuming it ran hot — and only calling what's left over RNA.

### 2. Bayesian posterior credible interval (principled)

In the full mixture model, the posterior distribution over n_RNA given the data is:

> P(n\_RNA | k, n) ∝ P(k | n, n\_RNA) × P(n\_RNA)

where P(k | n, n_RNA) is the likelihood of observing k=800 given that n_RNA molecules are RNA and n − n_RNA are DNA, integrated over both beta-binomial sampling processes.

You compute this posterior numerically — it's just a sum over all plausible values of n_RNA from 0 to 1000 — and report the **5th percentile** as your lower credible bound.

---

## The Uncertainty Has Two Sources

This is the key subtlety — the width of your credible interval comes from two places:

**Sampling noise**: Even if you knew exactly n_RNA, the beta-binomials produce stochastic counts. More overdispersion in the DNA component means wider uncertainty.

**Component ambiguity**: Some molecules genuinely can't be assigned — if the DNA and RNA beta-binomials overlap near, say, 60/40 splits, those molecules are inherently uncertain.

The credible interval properly propagates both.

---

## Practical Consequence

For an 800/200 example, the separation between μ_DNA=0.5 and μ_RNA=0.99 is large, and n=1000 is substantial. In that regime, the posterior on n_RNA will be fairly tight, and the 5th percentile won't be far below the mean estimate. But for a region with n=50 and a 35/15 split, the credible interval will be wide — correctly reflecting that you can't confidently deconvolute with sparse data.

This gives you a natural **per-region quality filter**: regions where the lower credible bound on n_RNA is near zero get flagged as "ambiguous," while regions where even the conservative bound is large get called confidently RNA-enriched.






### splicing anchor tolerance

fragment overlaps regions
when overhang into a region is tiny <=3bp, clip
avoids tiny overlaps due to misalignment
need to reduce total read length so mass is conserved



### Build and finalize FL models

- Generate mini-genomes
- Generate arbitrary transcripts over mini genomes
- Generate simulated fragments with different RNA and gDNA FL distributions. try multiple combinations of fragment lengths, including where the RNA FL is much greater than gDNA FL and vice versa. RNA FL >> gDNA FL and RNA FL << gDNA_FL. Try when they are equal.
   - Set RNA and DNA FL profile
   - ex. RNA FL mean 150, stdev 20. DNA FL mean 70, stdev 20
   - Simulate fragments over mini genome
   - Run rigel
   - Measure RNA and gDNA FL estimates, measure error
- Determine if FL estimation is accurate in simulation


### region density estimation

Each region has 'contained' and 'boundary crossing' fractional counts.

Estimating the 'contained' fragment density is simple:

density_contained = contained_fragments / (region_eff_len)

- contained_fragments: the fully contained fragments in the region
- region_eff_len: FL corrected region size. But which FL distribution do we use? we have mixture of gDNA and RNA.

gdna_density_contained = gdna_contained_fragments / (region_eff_len using the gDNA FL dist)

rna_density_contained = rna_contained_fragments / (region_eff_len using the RNA FL dist)

Next, I want to *model* the density rather than just compute it for each region.

Let's start with INTERGENIC regions. We need to build a model that can predict: 

number of fragments ~ region size

Or probably should model:

log(num fragments) ~ log(region size)

This is a density model right? Because density = num_fragments / region_size

But if we try to build a model of this across all regions, we can learn the pattern of density.

How do we model this? 

log(num fragments) ~ log(region size)

My thought is quantile regression so that we can then have a model that we can query:

Given a region of size L bp, what's the 95th percentile number of fragments I would expect to see?

What are other ways of modeling this? 

I want to build this model for INTERGENIC contained fragments and also for INTRONIC contained fragments (intron is mixture of nascent RNA and DNA). For EXONs which are a mixture of gDNA nascent RNA and mature RNA, we can then predict RNA density

- given exon region of size L bp, how many are gDNA (95th percentile)?
- subtract that amount
- whatever is left is RNA

Can you audit the current code. Help me design and plan this kind of model. This will work for unstranded data.

We still need to incorporate boundary crossing fragments into the model



For INTERGENIC regions, we use the gDNA FL distribution because we assume INTERGENIC ~ gDNA. For INTRONIC regions we assume are highly gDNA enriched we can also assume the gDNA FL distribution for now.

We need a generic effective length estimation system. Given the total region length (region_end - region_start) and a FL distribution, we should be able to estimate the effective region length. I want to make sure our methodology for this is extremely stable and robust. I believe it is but audit this and report your audit. Do we need to refactor or rename this? I believe we have this in place I just want it to be clean, concise, elegant code that is named appropriate, algorithmically correct, and easily interpretable.




This gives us an initial INTERGENIC and INTRONIC density estimate for fully-contained fragments. This can form our initial gDNA global density estimate.

However, using intergenic and intronic regions to estimate gDNA density *FAILS* when we have hybrid capture data, because fragments will be concentrated on exons. This is where we will utilize 1) strand-specific data and 2) boundary-crossing fragment data.

**We do need to implement mappability-corrected length for gDNA estimation, because a substantial fraction of the genome is not mappable and should not be used in the denominator for gDNA density estimation (gDNA density will be underestimated until we correct for mappability)** We have the opportunity to incorporate mappability data from the 'alignable' tool (a tool that I built) as part of the rigel index build. Are there any relics of this still present in the rigel index code? This was implemented at one time but might have been removed. This can be re-implemented in a future enhancement. For now, we will get the nuts and bolts of gDNA estimation working.

The first step for gDNA estimation is standard density estimation for intergenic and intronic regions. Any fragments in these regions are going to gDNA enriched and mature RNA-depleted.


### density modeling



### gDNA vs RNA 

When we iterate over regions we can compute the density of the contained fragments.



### region strand data

When we do gDNA estimation, we need to model the *variance* of fragment distribution on the positive and negative strand. Biology mandates that gDNA is double-stranded and our mean is absolutely cemented and fixed at 50/50 positive/negative stranded. However, in my investigation of real data, the variance around the mean was not binomially distributed. There was real overdispersion. Modeling the strand variance will be important and something we will need to incorporate downstream. If variance is very overdispersed, we could see large strand asymmetry (for example, 10 total fragments, 8 sense, 2 antisense) and still have it be reasonably probable to expect it from gDNA. This needs to make it into the EM somehow -- perhaps by increasing the gDNA bayesian prior pseudocount or as a  regularizer (M-step?). 



## Remove boundary_kind instrumentation

This is likely not going to be used


## gDNA -> RNA failures


## coarse and fine variable name refactoring

Regarding variable names, I discourage names that involve "coarse" and "fine" region terminology. Discussing coarse and fine regions is solely for the purpose of planning. Once we implement this, regions are just regions and don't explicitly need to be labeled "fine" regions.


## gDNA FL hist

We have implemented Phase 0, Phase 1, and Phase 2 of our fine region migratino plan. We are preparing for Phase 3 and have created a planning document (attached) for Phase 3. I'd like you to review the document and help me improve it. There are a lot of good things there. I'm not quite convinced that the 'boundary_kind' fields are needed, and not quite convinced about exactly which 'boundary_kind' fields we should try to retain. I do agree that we need to decide which fragments we will use for gDNA FL estimation. That is crucial. Intergenic contained (signature 0x0) are useful for gDNA FL. Pure intronic contained (0x4, 0x8, 0xC) should go to gDNA FL. Then, we need to include boundary-crossing fragments. Let's look at how to include boundary-crossing INTERGENIC and INTRONIC fragments with the new system. Given a boundary-crossing fragment, we would only want to accumulate the INTERGENIC/INTRONIC portion of the fragment. For example a 100bp fragment overlaps an exon (right boundary, 98bp) and an intron (left boundary, 2bp). This is an exon-intron boundary crossing fragment. When we accumulate the fragment weight 2/100 on the intronic LEFT boundary, we would add this as a weighted FL estimate to the INTRON gDNA FL estimation. The intronic gDNA FL estimate will be weighted by 0.02. I would say it would be reasonable to maintain the following FL histograms for gDNA FL estimation:

- INTERGENIC contained
- INTERGENIC boundary-crossing
- INTRONIC contained
- INTRONIC boundary-crossing

Granted we will likely aggregate these gDNA FL histograms into one, but it makes a lot diagnostic sense to keep them separated during accumulation to make sure the FL hists are similar/compatible and how much they are being contaminated.

We can maintain the other FL hists:
- EXON contained
- EXON boundary-crossing

Initially we will not use the exon FL hists for gDNA FL estimation. There is a possibility that if we partition regions into 'expressed' (RNA+) and 'not expressed' (RNA-), we could then utilize fragments in the 'not expressed' regions, which should be relatively pure DNA. This was a previous strategy. But we might not need to do this because we tend to get many fragments in the intergenic/intronic contained and boundary-crossing categories. This is a future endeavor.

So I think the next step is to refine the plan for gDNA FL estimation as above. Evaluate my design above. What's your critique? Is this is a good design. Improve the implementation and Phase 3 to gather FL histograms appropriately.

What are the implementation steps ahead of us?

What interfaces and design decisions do I need to make to utilize our new region partition and fractional accumulation in downstream code. What are the downstream consumers besides the FL model.



Re: subtleties to nail down 1) Yes! gDNA FL pools are unspliced fragments only. 2) Agree that gDNA FL pools collapse pos + neg strand (DNA is unstranded). 3) there are lots of 'ambiguous' signatures that involve exons + introns. these are impure regions that do not give us a clear signal. Now, we have the opportunity to separate out the ambiguity further. INTERGENIC is trivial (0x0). INTRON "pure" (intron_neg=0x4, intron_pos=0x8). INTRON "ambig" would be just intron_pos+intron_neg=0xC. This concept translates to exon regions too: EXON "pure"  (exon_neg = 0x1 and exon_pos=0x2) and EXON "ambig" (exon_pos+exon_neg=0x3). The rest of the signatures are "mixed" because they have overlapping EXON and INTRONs. While I may have suggested that we consolidate our categories, having thought this through, I think it is going to be helpful to 




## increasing region granularity

One of Rigel's main goals is to correctly deconvolute genomic DNA (gDNA, unspliced, unstranded), nascent RNA (nRNA, unspliced, strand-specific), and mature RNA (mRNA, spliced and unspliced, strand-specific).

We continue to have higher than acceptable error. This PR attempts to decrease pool-level error, correctly resolving gDNA, nRNA, and mRNA.

The tool partitions the genome into non-overlapping regions. The current behavior is to cluster overlapping exons into a single conglomerated region. This mostly works well, but occasionally yields some large exonic regions that contain multiple overlapping transcripts and genes. Sometimes the culprit is a very long single-exon transcript. Other times it is can be read-through transcripts that splice from one gene into another. Whatever the cause, large exons destroy granularity that we need. We hypothesize that restoring granular regions (instead of coarse regions) will improve the tool behavior overall.

This PR is to refactor our region partioning and the downstream analysis. This PR will touch many parts of the code including the Rigel index, BAM scanner/buffer, scoring, calibration, locus building, and setting the bayesian prior and warm start for EM.

Fortunately, Rigel USED to have a granular region system in place but we abandoned this in favor of the current coarse but simpler region system.

You can audit and retrieve code for the prior region system from the git repo as commit '14c307f'.

Our new granular region system, which I want you to reimplement, is defined as follows:

A region is defined as a genomic interval with (ref_id, start, end). We will partition the entire genome into non-overlapping regions.

Each region is defined by a set of true/false flags:
- intron_pos: an intron on the genomic positive strand overlaps this region
- intron_neg: an intron on the genomic negative strand overlaps this region
- exon_pos: an exon on the genomic pos strand overlaps this region
- exon_neg: an exon on the genomic neg strand overlaps this region

This 4-bit, 4-flag signature defines a region. Note: we may have used tx_pos and tx_neg instead of intron_pos, intron_neg in the past. These are interchangable; one can be computed from the other.

if bits are intron_pos, intron_neg, exon_pos, exon_neg

b0000 = intergenic
b0001 = exon neg
b0010 = exon pos
b0011 = overlapping exons on pos and neg strands
b0100 = intron pos
b0101 = intron on pos strand overlapping with an exon on the neg strand
b0110 = intron on pos strand overlapping exon on pos strand

.... and so on and so forth hopefully you understand the pattern here?
All combinations of 4 bits are possible


The region partitioning can then be defined. Given a genome (fasta file) and gene annotation (gtf file):

- For each reference chromosome:
   - create a sorted list of unique exon boundaries by iterating over transcript exons. each boundary has the following info: (is_tss = is the boundary of a transcript start site, is_tes = is the boundary of a transcript end site, is_exon_start, is_exon_end)
   - initialize 4-bit state as b0000 (intergenic) and iterate from left to right over sorted boundaries
      - if state *changes*, then save current region and its state; update 4-bit state for each boundary encountered

This is roughly how the region partition definition algorithm works.

Today's "coarse" region partitioning can be easily derived from the fine-grained region partitioning.

The coarse partioning is unstranded and only has INTERGENIC, INTRON, and EXON classes. We also have an ambiguous (AMBIG) class.

if bits are 3:intron_pos, 2:intron_neg, 1:exon_pos, 0:exon_neg

4-bit states new fine-grained system -> old coarse system
b0000 = INTERGENIC
b0001 = EXON
b0010 = EXON
b0011 = EXON (AMBIG)
b0100 = INTRON
b0101 = EXON
b0110 = EXON (AMBIG)
b0111 = EXON (AMBIG)
b1000 = INTRON
b1001 = EXON (AMBIG)
b1010 = EXON
b1011 = EXON (AMBIG)
b1100 = INTRON (AMBIG)
b1101 = EXON (AMBIG)
b1110 = EXON (AMBIG)
b1111 = EXON (AMBIG)

You'll need to verify and derive this for yourself. This is meant to give you a pattern. Just looking at this table should make it clear why the coarse-grained region partioning obfuscates a lot of important granularity. 

We need to revert our code base to the fine-grained region system. This will give us finer regions to evaluate gDNA density. It should yield better estimates of gDNA density. It should yield better estimates of exposure weight for gDNA effective length normalization.

I need you to do the following:

1) find code in the older git repo that contained the prior region calibration system
2) plan a migration from coarse-grained regions to fine-grained regions. This will need to happen in phases.

Here are proposed phases:

Phase 0: code cleanup. Audit the rigel index, region partioning, bam scanning, and calibration code. Try to clean up aspects that are stale, overly complex, legacy, dead, etc. We want to implement over a clean slate.

**NOTE** We are going to implement this with ZERO backwards compatibility and ZERO legacy compatibility

Phase 1: update rigel index to use fine grained regions. rewrite the rigel index code to derive fine-grained regions instead of coarse-grained regions.

Phase 2: BAM scanning. In the BAM scanner, we accumulate fragments over matching regions. We are going to revamp the way we do this accumulation. The new accumulation will give us powerful new information that will be used during calibration.

When we are scanning the BAM file, we construct Fragment objects. A fragment object is a collection of aligned blocks (ref_id, start, end) and other information. We need to map fragments onto regions. A lot of the code to do this can be reused.

Currently, we exclude multi-mapping, chimeric, splice artifacts, and perhaps other classes of transcripts that are difficult to interpret during the calibration step. Likely we can use the existing fragment exclusion criteria.

For region accumulation, each regions needs to store fragment count data as fractional counts (float). Fragments are either CONTAINED (completely contained within a single region) or CROSSING (overlapping multiple regions, e.g. generates boundary flux)

- CONTAINED unspliced pos strand, unspliced neg strand, spliced pos strand, spliced neg strand
- BOUNDARY FLUX LEFT unspliced pos, unspliced neg, spliced pos, spliced neg
- BOUNDARY FLUX RIGHT unspliced pos, unspliced neg, spliced pos, spliced neg


For eligible fragments, we map them over the region partioning using the compatible aligned (ref_id, start, end) blocks.

- For each fragment (spliced and unspliced):
   - for each aligned block (ref_id, start, end):
      - determine the set of overlapping regions
      - if overlaps a single region update region data for contained
      - else for each overlapping region R_i (region_ref_id, region_start, region_end):
         - compute number of overlapping bases (overlap_bp)
         - fractional overlap is overlap_bp / aligned_block_size (end - start)
         - determine whether left, right, or both boundaries are crossed
         - update boundary flux for each boundary crossed. fractional overlap must be divided by number of boundaries crossed (no double counting)

We required that each fragment produces exactly 1.0 count units that are partioned among contained and boundary crossing compartments. For example, take regions R1 =(0, 100), R2=(100, 150), R3=(150, 300).

Fragment has aligned block (80, 140). Aligned block is 60bp Fragment overlaps R1 and updates R1 boundary right with (20/60)=0.33, then R2 boundary left with (40/60) = 0.66. 

Another example fragment with aligned block (50,250) overlaps R1 right boundary (50/200) = 0.25, R2 both left and right boundaries so R2 left (50/200/2) = 0.125 and R2 right = (50/200/2) = 0.125, then R3 left boundary (100/200) = 0.5. The sum of fragment contributions = 1.0.

This accumulation paradigm accomplishes the following:
- fine-grained
- maintains strand information
- maintains spliced/unspliced information
- allows computation of boundary flux (compatible with today's behavior)
- allows regions to compute their "independent" fragment support by aggregation boundary crossing and contained fragments

I believe that this region accumulation system will support a more sophisticated calibration step.

We can recover existing behavior from the region data if we choose to. But we can also extend existing behavior. We need to develop enhanced calibration that uses fine-grained regions and boundary flux.

Phase 3: Fragment length model

During BAM scanning, we must continue to identify gDNA-compatible fragments for the gDNA FL model. We can keep this approach the same as the current approach (INTERGENIC, INTRONIC, and EXON-INTRON fragments are used for the gDNA model).

We may choose to revert to an older gDNA FL model when we get here.

Phase 4: Calibration

This is not designed yet.


Begin designing this new implementation. Store the design in 'docs/fineregions'. 



## argument for fractional accumulator

I want to discuss fractional 12-channel accumulator that you felt that we should defer indefinitely. Looking through your rationale now -- 1) It is wonderful that the current CalibrationAccumulator is sophisticated. We would KEEP most of the current infrastructure when we change to float32 (not float64) fractional counts. We would expand boundary flux with orientation to include spliced fragments.  2) We would KEEP splicing anchor tolerance (excluding fragments from flux when overlap is less than `tolerance` in bp). The current splicing_anchor_tolerance support is wonderful and should work. 3) Integer accumulator creates a "double counting" problem -- suppose a fragment overlaps 7 fine-grained regions. It would increment boundary flux in all regions, but there would be no way to know that the incremented boundary flux came from a single fragment. We have fragments of many different lengths: an 800bp fragment might overlap many more regions than a 150bp fragment. This creates bias in integer accumulation to longer fragments. For coarse-grained accumulation, we did not pay attention to this but with fine-grained accumulation, we are breaking apart larger EXON regions into many smaller regions. Fragments will often overlap many regions. Why do you feel so strongly about deferring our proposed fractional 12-channel accumulator? This concept would allow each region to independently accumulate fragment support. One of the nice properties is that for a single region, we can sum the boundary flux and contained fractional support in a region and easily obtain a total. The sum across regions will equal the total fragments used in accumulation. Without this, we cannot easily solve "how many fragments are in this region?". Right? ---> In summary, I don't see or understand how an expanded fractional accumulation creates major problems. This would add support for spliced fragment boundary flux (new feature) AND make boundary flux universally helpful rather than just used for EXON-INTRON boundaries. We would use boundary flux everywhere, for all regions. -----> I have not developed the algorithms yet, but I believe that this expanded fractional accumulator formulation will pave the way for more sophisticated gDNA estimation. I am thinking about complex sets of regions with overlapping introns and exons and both pos/neg strands. Instead of one coarse EXON region, we will have many adjacent exonic regions. My plan is to devise a gDNA density estimator for these internal exon regions that incorporates spliced/unspliced boundary flux and contained fragments to predict gDNA densities. I don't think this would be possible without this formulation. If you have a better formulation or can suggest improvements please do! Also, I would love you to do make the implementation extremely efficient! We could certainly boost performance if we think about implementation ideas. This is a first pass at giving us the foundation we need.

### splicing anchor tolerance goes away

Splicing anchor tolerance logic needs to be examined in the setting of our new fractional accumulator, because it may no longer be needed. The splicing anchor tolerance was intended for fragments that cross exon-intron boundaries (exon region -> intron region) but the intronic overlap was very small (<3 bp default). this small "overhang" could be aligner error where that last few bases may truly belong on the next exon (spliced) but the aligner didn't align it correctly. when we are computing exon-intron boundary flux, we want to use that as a proxy for gDNA and the fragments with tiny overhang into the intron would otherwise contaminate our exon-intron signal with RNA. However, let's rethink this because the new fractional accumulation paradigm may elegantly solve the splicing overhang problem without the need for a tolerance. In the new fractional accumulator, a fragment that overlaps an exon-intron region pair and overhangs into an intron would contribute to boundary flux, but that flux would get partitioned between the exon and intron. If the aligned block is 100bp with 98bp exon and 2bp intron overhang, the flux would be partitioned so that 98/100 would be allocated to the exon and 2/100 would be allocated to the intron. So if we use the intron region boundary flux to compute gdna density, we would only see 2/100 contribution from the overhanging read -- it is true that this is rna noise, but far less than the prior case where every boundary crossing fragment gets equal integer weight of 1. The intron-to-exon flux which will be useful for gdna estimation is weighted by intronic overlap. So it may be the case the the splicing anchor tolerance isssue is naturally handled by this approach. Can you evaluate this carefully? Derive the logic yourself and work through some examples. Can we rely upon our fractional accumulator to handle splicing anchor tolerance naturally? If so, this simplifies our algorithm. Go ahead and rework the design document v3 to remove splicing anchor tolerance logic, unless you believe that it still serves an important purpose

## hybrid capture simulator

We have an extensive code base for simulation in place in the 'src/sim' folder within Rigel. To support hybrid capture simulation, the simulator needs additional inputs.

1) The inputs are 'probes' which can be provided in transcript coordinates (transcript_id, start, end) or genomic coordinates (BED12 format in case probes span multiple exons and involve multiple genomic intervals). 

2) Hybrid capture changes the probability landscape for sampling reads. The current model samples uniformly along a transcript (or the genome). A capture model changes this to a non-uniform probability landscape.

We need a simple and straightforward probe-binding hybridization energy model to create a probability distribution. The energy is related to the number of overlapping (hybridized) bases along the probe. Can you research how this model should work? A simple idea linear in the number of overlapping bases between a fragment and the probe. 

Here's a worked example:

Given transcript T with length 2000bp. It has 3 probes at (200,320), (500,620), and (1500,1620).

- A 200bp fragment at (1000,1200) overlaps no probes and gets no hybridization energy. It is effectively off target
- A 200bp fragment at (1310,1510) overlaps a probe by 10bp and gets +10 energy. 
- A 200bp framgent at (1480,1680) completely overlaps the probe and gets full +120 energy.

Can you research how to implement this? Ultimately, we need to create a probability landscape along a transcript that we can sample from. Sampling from regions that overlap probes will be much more likely than sampling from regions that do not overlap probes. Probability needs to be relative to the overall transcript abundance (relative to other transcripts).

Transcript probes "project" onto the genome. A transcript probe may cross an exon-exon junction and be "spliced" such that its genomic project is split into multiple genomic intervals.

When we sample gDNA (for simulating gDNA), the gDNA sampling must use a similar energy model, but when probes are split, the potential binding energy is much lower, so the RNA transcripts have a big advantage compared to when the transcript probe is on a single genomic span.

We need the hybrid capture simulator to have the following properties:
1) simple (not too complicated)
2) easily implemented
3) roughly biologically accurate (doesn't need to be perfect.. just "good enough")
4) extremely good performance -- we don't want to be waiting hours to simulate reads. this needs to be blazingly fast.




## gdna capture weighting

I'm just thinking out loud here.

We already have regions partitioned.
- Intergenic (gdna)
- Intronic (gdna + nrna)
- Exon-intron (gdna + nrna)
- Exon (gdna + nrna + mrna)

Given that nascent RNA tends to be very sparse genome-wide, intergenic + intronic + exon-intron are reasonable to approximate gDNA density.

Each region is a genomic interval and after this step we have gDNA density estimates for each region. EXON regions are also estimated by projection of the bordering exon-intron regions. So all regions should have gDNA estimates. That should already be happening in the code.

The simplest next step is to figure out how to weight each region. Could each region could be assigned a weight?

What if we divide/partition the total gDNA "mass" among the regions? Each region gets a weight proportionate to its contribution to the total gDNA mass.

So for example let's say there is the 100kb locus with 1kb exon and 99kb intron. We have the exon gdna densities (computed from the exon-intron boundaries) and intron gdna densities. There are 1000  total gDNA fragments in the 100kb locus, but 950 of them are consolidated on the exons (950/1kb = 0.95 frags/bp density) and 50 of them are in the introns (50/99kb = 0.0005555 frags/bp density).

We then would assign per-region weights such that the GLOBAL sum of the weights equals 1.0, and the regions are then weighted by their contribution to the total global gDNA mass.

This seems extremely simple to me. The aggregate gDNA locus level estimates are then multiplied by their weights so that the high-density regions contribute much more weight than the low-density regions.

I do recognize that there is a "circular" logic issue here where we are weighting gDNA estimates by their densities.. kind of but not exactly effectively using density squared as the weight rather than density alone. On that note, raising the density to a power effectively squashes small values and retains larger values, effectively distributing density across more highly-weighted regions.

You mentioned an idea to subsample the total region population to build a density distribution so that we are using a portion of the total regions to construct a density distribution and then using that distribution predict per-region weights. I'm not sure if that helps but repeated subsampling and computing weights would incorporate some aspect of uncertainty in the region measurements.

Speaking of uncertainty, we already had a plan to enhance how we compute bayesian gdna prior. The current bayesian gdna prior does not include uncertainty in any way. It is computed as a flat pseudocount that adds. The gDNA estimation has an uncertainty component that ideally should be modeled.

Uncertainty comes from several sources. First, tiny regions have much greater uncertainty than large regions (big difference in exposure). Exon-intron boundary regions are small and there will be substantial uncertainty there.

Uncertainty comes from not-modeled factors that create gDNA variation such as a GC content.

When we apply gDNA weighting, we should consider incorporating uncertainty modeling at the same time.







We don't have enough information yet to be able to change our code or parameters. We need to delve deeper into this issue.  FIRST, I think we can isolate the problem by eliminating/discarding multi-mapping reads. First audit the code, there should be a CLI parameter include multimapping that turns on/off inclusion of multi-mapping reads. Interrogate the code flow and ensure that this is implemented correctly. IF we turn OFF multimapping reads, we should break the mega-locus issue. If the mega-locus problem is the main contributor to gDNA -> RNA errors, we should be able to diagnose and visualize this after discarding multimapping reads. Without multimapping reads, we should see a huge performance boost. Please do this first. 

## ambiguous region errors

We don't have enough information yet to be able to change our code or parameters. We need to delve deeper into this issue. I want you to find regions where gDNA fragments are mispredicted to be RNA using the annotated BAM output file. Select up to 5 regions that are particularly egregious. Show that these regions are so we can look at the BAM file on the genome browser and a get a sense of what the reads look like -- find precise locations where there are a very high number of mispredicted fragments. Analyze these on an individual per-fragment level to understand what is happening to the likelihoods of the fragments. We need to find out exactly what is happening to these fragments. 


## mappability corrected effective length
   I agree strongly with your concern. Intergenic calibration without mappability-corrected effective length is fragile, and hybrid capture breaks naive intergenic extrapolation. The denominator should become something like accessible/capturable effective length, not raw genomic span. For capture libraries, gDNA exposure should include probe/bait targetability or an empirically learned capture profile.

## gdna prior uncertainty
  **Treat the gDNA prior as uncertain, not as a hard pseudocount.**
   Right now the calibration-derived `gdna_prior_count` enters as a physical Dirichlet count. In unstranded data, especially when calibration channels may contain nascent RNA, that prior should be downweighted by an evidence-quality factor. Good knobs to test: cap per-locus `eta_g / n_em_fragments`, shrink the prior toward zero when strand contrast is unidentifiable, and propagate calibration uncertainty instead of using a point estimate.

## gdna vs nrna eff length normalization
   Current gDNA length normalization is not totally unconstrained; it is baked into the gDNA log likelihood using a local sampling window based on transcript span plus a gDNA flank. Synthetic nRNA then gets transcript-like effective-length normalization in EM. For long gene spans, the nRNA “shorter than gDNA” advantage is probably tiny, and for exon-contained fragments mRNA has the real effective-length advantage. I would audit whether gDNA should use the full locus/capture-accessible union instead of an anchor transcript span, but avoid arbitrary penalties that create false nRNA.


## Remove 'coverage_weights'?

Not used as prior anymore
Used as warm start
But does it do anything?


## Rigel index with GTF duplicates

GENCODE GTF has duplicates (why?)
Need option to drop duplicates during index build
Keep lexicographically smallest transcript id?

## Performance optimization


## Locoregional calibration

Currently, calibration collect global information to estimate the gDNA prior for bayesian EM.  This might be improved by incorporating locoregional information, which may reflect copy number changes.

Given a particular Locus (genomic region containing transcripts) this approach would compute gDNA densities from the flanking regions of some genomic distance (how much? 10kb, 100kb? 1Mb?). The locoregional gDNA estimates would shrink to the global estimates using empirical bayes shrinkage and then feed into the bayesian prior for the Locus of interest.



# Polars migration

Polars is apparently much faster than pandas. Might be ideal for rigel

## cgranges interval searching

### possible to serialize cgranges index? saves time and memory...

### should we "lazily" build cgranges index if only used during one stage of pipeline? for example, calibration requires cgranges index for mappability but then might not be used again. could lazily build and free, potentially decreasing overall RSS







## Gene counting

Ensure gene quantification is correct.

Regarding gene quantification. Yes, we can sum counts of the transcript isoforms of the gene. Annotated single-exon transcripts are associated with a gene and will be included in those counts (they are both mature and nascent so this is okay). For gene quantification, we MUST exclude synthetic nascent RNA (synthetic nascent RNA should not be associated with a gene. Remember, we can sum counts, but when we compute transcripts-per-milion (TPM), we use a weighted average of counts respecting different effective lengths.




## Transcript versus pseudogene multimapping (edit distance improvement)

It appears that transcript vs processed pseudogene identifiability is the largest remaining source of error in rigel. When running with minimap2, rigel fails to distinguish the true transcript from its processed pseudogene, resulting in massive counting error.

Our hope is that the 'true' transcript will contain fewer mismatches/edits and can be identified by superior alignment quality compared to other competing alignments.

We wish to improve how we distinguish multimapping alignments. We want to be able to prioritize the 'true' alignment relative to 'false' alignments.

The rigel tool currently tries to do this by parsing the 'NM' tag in the BAM file. This tag is supposed to represent the 'edit distance' between the read and the reference. Using 'NM' is intended to help distinguish alignments by the number of edits between the read and the reference. However, some aligners may or may not be setting the 'NM' tag appropriately. Some tools only set 'NM' when there are insertions/deletions. It is unclear how minimap2 sets the 'NM' tag.

The 'MD' tag in SAM/BAM files is a companion to the CIGAR string and is intended to identify the exact positions of the mismatches/edits. 

Here's an excerpt from the PDF describing SAM/BAM optional fields, that defines and explains the MD tag:

### MD tag definition

MD:Z:[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*

String encoding mismatched and deleted reference bases, used in conjunction with the CIGAR and SEQ fields to reconstruct the bases of the reference sequence interval to which the alignment has been
mapped. This can enable variant calling without requiring access to the entire original reference. The MD string consists of the following items, concatenated without additional delimiter characters:
- [0-9]+, indicating a run of reference bases that are identical to the corresponding SEQ bases;
- [A-Z], identifying a single reference base that differs from the SEQ base aligned at that position;
- \^[A-Z]+, identifying a run of reference bases that have been deleted in the alignment.

As shown in the complete regular expression above, numbers alternate with the other items. Thus if two mismatches or deletions are adjacent without a run of identical bases between them, a ‘0’ (indicating a 0-length run) must be used to separate them in the MD string. Clipping, padding, reference skips, and insertions (‘H’, ‘S’, ‘P’, ‘N’, and ‘I’ CIGAR operations) are not represented in the MD string. When reconstructing the reference sequence, inserted and soft-clipped SEQ bases are omitted as determined by tracking ‘I’ and ‘S’ operations in the CIGAR string. (If the CIGAR string contains ‘N’ operations, then the corresponding skipped parts of the reference sequence cannot be reconstructed.)
For example, a string ‘10A5^AC6’ means from the leftmost reference base in the alignment, there are 10 matches followed by an A on the reference which is different from the aligned read base; the next 5
reference bases are matches followed by a 2bp deletion from the reference; the deleted sequence is AC; the last 6 bases are matches.

### Double-counting edit distance when paired-end reads overlap

- When a fragment size is less than 2X read length, then the paired-end reads sequence a common portion of the fragment. For example if the fragment length is 250bp and we do 2 x 150 bp paired-end sequencing, then the two reads redundantly sequence a common 50bp region of the fragment. When we calculate edit distance, we have to ensure that mismatches/edits in this common "double sequenced" region do not get double counted.

- The STAR aligner avoids this double counting with the 'nM' tag
- Minimap2 and other aligners do not appear to handle this issue


### Proposed solution and implementation steps

- Instead of using the annotated BAM file, let's start by creating a simulated scenario using the GAPDH gene. 
- We have a wonderful 'sim.py' script in 'scripts'. We can set this up to produce reads to GAPDH transcripts (and zero transcripts everywhere else). Then we can generate a tiny simulated scenario that can be used for debugging.
- The simulation needs the full human genome reference so that reads multimap to other locations.
- Generate simulated paired-end reads to transcripts from the GAPDH gene. GAPDH has many pseudogenes and reads often multimap. The current rigel tool with minimap2 has massive errors related to GAPDH abundance estimation.
- Create a simulated "oracle" BAM file to GAPDH
- Align the reads using minimap2
- Run rigel on the oracle BAM 
- Run rigel on the minimap2 BAM
- Trace the reads/fragments through both the oracle BAM and the minimap2 BAM.
- Find fragments that multimap. Understand how these are being processed. Understand how the likelihoods are calculated. Understand at a deep level why Rigel cannot correctly distinguish these reads.

#### Big question: can we distinguish transcript from pseudogene using any aspect of the alignments?

- NM tag?
- MD tag?
- Other alignment results?
- How would we design this?



## Genome to transcript coordinate search

This code may be duplicated/redundant in more than one place? It is being used in scoring and also in fragment length calculation. Can we consolidate code? Can we improve the efficiency of data structures?


## Unannotated splice junctions

A fraction of splice junctions are "unannotated" in that they don't have exact matches to the reference. These can be either artifacts (aligner errors) or 'novel' isoforms/transcripts. There is a lot of information that we are not taking advantage of yet:

- partially annotated: 3' site, 5' site, or both
- (easy) is the 5' splice site annotated
- (easy) is the 3' splice site annotated
- (medium) how much 'anchor' on either side of the intron?
- (hard) does this junction match any of the 'blacklist' splice junctions -- requires building a splice junction blacklist.






## Equivalence class sorting?

This seems to affect non-deterministic behavior of the EM. Is this necessary?

    // ---- Deterministic ordering ----
    // Multi-threaded BAM scanning produces fragments in non-deterministic
    // order.  The unordered_map above inherits that non-determinism in both
    // (a) the iteration order of equiv classes, and (b) the row order of
    // units within each class.  Since the EM E-step accumulates column sums
    // over rows, and FP addition is non-associative, different row orders
    // produce ULP-level differences that SQUAREM amplifies across iterations
    // potentially causing large cascading output differences.
    //
    // Fix: sort equiv classes by comp_idx, and sort rows within each class
    // by their log-likelihood fingerprint.  This makes the EM iteration
    // fully deterministic regardless of input fragment order.

