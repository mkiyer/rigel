Model A -- 

At initialization we do not know what the gDNA landscape will look like (hybrid capture on/off? Cancer sample with lots of amplifications/deletions?). We must have an unsupervised approach.

Unsupervised means we can only start with observable nodes. At initialization, these are just intergenic and intronic nodes. I would recommend a log-log LOESS (Exactly what you described). This needs to be built from intergenic regions, intronic regions, exon-intron boundaries, and exon-intergenic boundaries. Please create and plot this so we can see how it looks and set the span appropriately. If we don't have enough data (regions and boundaries) for a proper LOESS fit, we should create a larger simulation dataset. My preference is to try to keep the existing dataset to facilitate rapid iteration and prototyping (having to run many large scenarios slows development time). 

I agree that our variance ~ mean models have become quite disorganized and messy. Perhaps we can consolidate them. Here is my understanding.

1) DIRECT observation variance ~ mean model (this model). This models the variance ~ mean relationship directly observed at nodes. Expect a bimodal pattern under simulated hybrid capture conditions. For real rna-seq data, we expect greater dispersion/variance should be greater and less predictable. Our variance ~ mean fit should be monotonic. We have LOESS in mind but this could be substituted for a different variance ~ mean fit that can respect the bimodal or otherwise nonparametric behavior that we see (allows the data to dictate the fit rather than imposing a rigid distribution). I'm flexible if other models could work better.

2) IMPUTATION variance - This model is built from adjacent boundary pairs (or boundary-region-boundary triplets). This currently is a var ~ mean log-log LOESS. Similar constraints - should be monotonic. Should learn the data pattern and compensate for enrichment (hybrid capture) and depletion, when present.

A challenge with these models is that hybrid capture panels sometimes enrich thousands of genes (whole transcriptome) and sometimes enrich a very small number of genes (targeted panels). Our approach needs enough locality to fit small targeted panels. A large span/bandwidth would "oversmooth" the fit when small targeted panels are present, losing the correct variance ~ mean relationship. Calibration of the locality of the fit is an open task. A single span parameter is not appropriate for all datasets. The span must be calibrated to 1) monotonicity and 2) adapt to locality. If LOESS is not ideal, we could embrace a different method. We want the method to be simple and elegant. If in doubt, this could be something to research.

The IMPUTATION model could be built in different ways.

The region/boundary graph has a single path with region -> boundary -> region -> boundary.

A more generalized version of the var ~ mean fit could be:
- Iterate over node adjacent boundary pairs, forming (boundary <-> region <-> boundary triplets)

Then we have different cases:

- Boundaries observable but middle region is non-observable -> impute the region gdna from each of the two boundaries -> these two observations constitute different measurements of the same value. We just don't have the 'correct' answer.

- All three (both boundaries and region) observable -> measure the ability of the boundaries to predict (impute) the region. Here we have the region as a 'correct answer' and can model how well the boundaries can predict it.

The concept of an iterative solution is quite important. After the first pass through the tool, we will have estimates of all regions and boundaries. When we pass through a second time, we can build these imputation models. The caveat is that on subsequent passes the models will be built from imputation which may have imputation noise and from strand deconvolution which will introduce some noise. However, it will give us a large set of observations (much better than the seeds).

## Building a robust initial model

Building a robust initial model is challenging when we have hybrid capture conditions. This is for the reasons that you exactly highlighted above. Under hybrid capture conditions, our enriched regions tend to be in exons, and we don't have observable triplets with two boundaries and one region from which to train a imputation correctness model. The examples we do have are going to be most likely from depleted regions. So our initial model is going to be unlikely to correctly predict the imputation variance. It would be dangerous because the initial model won't have observed mean counts in the range of enriched regions if the observations come from depleted regions. and the initial fit won't have observed variance at the at the level of enriched regions.

There are a couple of options. We do have the opportunity to observe pairs of boundaries that anchor enriched regions when we have a intron to exon to intron pattern. So a pair of boundaries that measure a single region where the two boundaries are both observable. So we can model variance from pairs of observable boundaries, but we cannot have access to the triplet of data that includes the region itself that we're imputing. So we cannot really train a model to fit to the mean because we don't have a measure of the reliability unless we have access to both boundaries and regions.

The next idea is twofold. We can apply our solver to any trivial nodes that can be solved directly. In other words, if we have strand information or strand specific data, we could, as an initial pass, solve all of the trivial nodes that can be solved with strand specific data alone. that initial pass would give us access to boundary region boundary triplets and allow us to train a imputation model.

The extension in the second half of this thought is related to iteration. We already discussed that iterative solving would give the model fit access to previously nonobservable nodes. This is very similar to the idea of applying some form of solver and retraining the variance mean model fit. The key here is that our initial variance model in the absence of running our solver on any nodes, must be configured to be highly uncertain. In other words, when we first run this as a first pass and we have a lot of uncertainty, we cannot allow our initial models to give us a false sense of certainty about imputation and about count variance just because they were trained on reliable stable depleted seed regions. I think that is the danger.

Perhaps, we need to commit fully to the idea of an iterative model. This is the only way I can think of right now to resolve the circularity. Before we start to impute, we need a trustworthy imputation model. To get a trustworthy imputation model, we need examples to train with, to get examples to train with, we need to apply our solver to At least some of the nodes so that we have enough, uh, observations of the mean and variance behavior to to fit a model. Or we need to have Our first pass through the data happen with a very high variance first model. so that our initial estimations are appropriately uncertain and then can be refined in subsequent iterations.

In summary, we have the building blocks in place. The important building blocks are the ability to model the variance mean relationship, The Simplex Solver itself is a very important building block. Our region boundary boundary partitioning structure, all of the data accumulation, it's all ready. And now we need to find a way to apply it gracefully. In this turn, I want you to work on a way to do that, incorporating my ideas here.


# Initialization

First of all, I want you to write up a very detailed implementation plan document. We have "CALIBRATION_PLAN.md" but I fear that it is now out of date. We may need a fresh plan.

Second, and more importantly, I think I have a key realization that will help us make this method much more elegant. Here is my idea. In our initialization pass, we will allow our variance mean model fit to train on all nodes. In other words, we can iterate over sets of boundaries and regions. and we can initially treat the raw data before deconvolution happens as our initialization. This means that Our initial fit will take two boundaries and always try to use them to predict the count in the region. Obviously, without any deconvolution or solving, that initial region count is just the total count. The total count includes all unspliced fragments. So we would try to build a model that uses the two boundaries to predict the total region count even though we know that many of the regions will have lots of RNA. and we haven't done any deconvolution yet. So the initialization will be wildly inaccurate because the boundaries will not be able to predict the region's count very well at all when RNA is present. So the variance mean model should give us a very uncertain result at a particular mean. We should see very high variance.

I think this idea exposes a few interesting properties that are potentially exciting. The first is that we might now have a measurement of convergence or improvement from iteration to iteration. We're using two boundaries to predict a region. And as we iterate, we expect that that boundary data will more accurately predict the region. We have a measurement of prediction error now. and we can minimize that. That will give us insight into whether we have converged potentially. it may give us more insight than that, but I need to think about it more deeply.

Let's reframe this:

- At initialization, we make the assumption: "all fragments arise from gDNA".
- We fit our var ~ mean models on all nodes using the total unspliced fragment mass in each node. 
- We can measure the goodness of fit in various ways. Individual nodes can each have an 'error' (predicted variance vs actual variance as a ratio or difference)
- Allows us to see how well model is predicting individual nodes

Iterations:
- var ~ mean (will be initially horrible)
- solve
- iterate


Allows us to track convergence:
- Overall var ~ mean goodness of fit *should* improve over iterations
- Individual nodes should improve.


In terms of next steps, I mostly agree with your plan:
1) control flow should now be: fit var~mean -> solve -> pass-k refine -> converge
2) 
- DIRECT var ~ mean from node's *own* counts (ALL nodes, initially gDNA counts/mass == total unspliced counts/mass)
- IMPUTATION var ~ mean from boundary <-> region <-> boundary triplets (ALL nodes, initially gDNA counts/mass == total unspliced counts/mass)
3) The concept of 'confidence' for each node should be the observed variance minus the expected variance (variance predicted by the model) -- do you have other insights into how we can measure accuracy per node?

In summary, this new control flow resolves a few issues:

- No need for a special "Pass 0" with strand-only. Every pass is the same. We apply over solver as best we can using the available data (with or without strand specific data)
- Individual node confidence can be estimated from the fit error
- Stopping rule should now be governed by quality of our var ~ mean model fit

I would love your critique on this.

I would love if you could create a new design and implementation document that encompasses all of the changes we have recently made.


# Variance model refinements

## Model 1: global gDNA var ~ mean model (locus level observations)

This model is MISSING from our implementation.

- Treat each Locus as an observation. Assume gDNA is uniform within a Locus. Use the Locus nodes (regions and boundaries) to compute mean and variance for the entire Locus.

- How do we separate the problem into 'Locus' observations?

We have excluded multimapping reads, and so Locus can be defined by INTERGENIC regions (region signature defines this). We can partition a reference chromosome trivially by partitioning at intergenic region boundaries.

Note, that each locus should INCLUDE its anchoring boundaries. The anchoring intergenic <-> exon boundaries should be included.

- Intergenic *regions* can be grouped into another special observation (with or without intronic regions as well)

Currently we designate intergenic, intronic regions as 'seeds'. This largely goes away. 

However, by grouping by Locus, intergenic regions are left out (they are not part of any locus).

The intergenic regions themselves can be grouped into a single intergenic observation with mean and variance computed from all intergenic REGIONS (contained fragments).

If we wish to maintain a separate gDNA "floor" model, it can be modeled by combining all intergenic (contained) fragments and intronic (contained) fragments. This assumes the intronic REGIONS has minimal nascent RNA overall (nascent RNA is rare/SPARSE assumption)

### Procedure for model fitting

- Fit a gDNA var ~ mean model by combining observations across every Locus (intergenic regions comprise an additional locus)
- Use the same premise everywhere. At initialization, all unspliced fragments are gDNA. There is no RNA.
- the var ~ mean model fit may be quite poor for loci with abundant RNA


## Model 2: gDNA imputation variance model

We have for the most part built this model, but probably have not built it correctly. It will take more work.

The model iterates over adjacent boundary <-> region <-> boundary triples. For each triple, the 'region' count is considered "correct". Its two boundaries are considered "estimates".

We have two boundaries BL and BR which bookend region R. 

Constraints:
- At least one boundary must have unspliced mass > 0
- The region must have unspliced mass > 0

The DIFFERENCE between each boundary and region is a measure of imputation ERROR.

The procedure is:
- compute boundary gDNA density
- impute boundary gDNA density onto region
- compute difference between observed (region gDNA count) and estimated (boundary gDNA imputation). This is measure of imputation error.

We can MODEL imputation accuracy by fitting
region gdna density ~ boundary gdna density


## Model 3: RNA imputation variance model

Similarly, we have started building this model, but probably have not built it correctly.

The RNA imputation model is very similar to the DNA imputation model, except that the RNA model can only be fitted on specific boundary <-> region pairs.

Compatible boundary/region pairs:
- Be must at SPLICE JUNCTIONS
- The boundary must have spliced mass > 0
- The region should have unspliced mass > 0

Given a boundary B and region R, the DIFFERENCE between boundary RNA density and region RNA density is a measure of imputation ERROR.

The procedure is:
- compute boundary RNA density
- impute boundary RNA density onto region
- compute difference between observed (region RNA count) and estimated (boundary RNA imputation). This is measure of imputation error.

We can MODEL imputation accuracy by fitting 
region RNA density ~ boundary RNA density


## Minimization of imputation error

We can compute overall imputation error.

For each node, the total imputation error is the SUM of the gDNA imputation error and the RNA imputation error.

Iterating should minimize the total imputation error


## RNA is NOT the inverse of gDNA

We have independent imputation error models.

Every node is stores as {f+, f-, fg}

The RNA model imputes region RNA (on a single strand) from the stranded spliced fragment mass at the boundaries. The spliced fragment mass are independent fragments that provide independent imputation estimates.

For example, given transcript on positive strand with region R bounded by BL and BL with U unspliced fragments.

- At initialization, R === { f+=0, f-=0, fg=1.0 } because all 300 unspliced fragments are assigned to gDNA
- BL has its own { f+, f-, fg }.
- We impute BL spliced fragments `f+` onto R.
- Now we have `f+ imputed` from BL -> convert to counts
- We have `f+ current = 0.0` on R -> convert to counts
- Now we have RNA imputation error measurement
- Initially error will be huge because RNA+ and RNA- are set to ZERO.
- After SOLVING, R will have new fractions {f+', f-', fg'}


The total imputation error of a boundary-region pair requires knowing both the RNA imputation error and the DNA imputation error relative to the current region counts (RNA+, RNA-, gDNA).

We wish to minimize the total imputation error.


## Generalizing to node landscape and propagation

During propagation, we propagate linearly from boundary -> region -> boundary -> region -> etc.

If we start with a boundary <-> region pair.

Boundary B { f+, f-, fg }
Region R { f+, f-, fg }

Our goal will be to determine how to propagate between B and R. Signal propagation occurs left to right and then right to left, so we will propagate in both directions: Forward Boundary -> Region and then Reverse Boundary <- Region.

For an imputation B -> R, what information do we need?

Above, we described the error model. We can compute the DIFFERENCE between the B imputation estimate and the current R value.

Now, the problem becomes, how do we UPDATE R given imputation from B?

You are correct. Each node needs to store a measurement of 'confidence'.

Getting 'confidence' right is the CRUX of the entire problem. 


## Node confidence measure

Every node must store a current value: {f+, f-, fg}

Every node needs a measure confidence in its current value

