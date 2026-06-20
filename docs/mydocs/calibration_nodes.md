
# Proposed calibration workflow

## Overview of calibration

There are two types of nodes: Regions and Boundaries

Global inputs:
- gDNA FL dist
- RNA FL dist
- strand specificity
- RNA var ~ mean count imputation model
- gDNA var ~ mean count imputation model
- gDNA global density model (hyperprior)

Region Node inputs:
- signature (ex+, ex-, in+, in-) 4-bit flag
- contained counts (spliced+, spliced-, unspliced+, unspliced-)

Boundary Node inputs:
- L/R unspliced counts, L/R unspliced mass
- L/R spliced counts, L/R spliced mass

Imputation inputs:
- signal propagation from left -> right (fractions, uncertainties)
- signal propagation from right -> left (fractions, uncertainties)

Outputs:
- Each node must solve for its { f+, f-, fg } with uncertainties { unc+, unc-, unc_gdna }

Final outputs of calibration:
- Node fractions converted back to count space.
- gDNA counts feed the EM gDNA prior
- gDNA feeds the effective length calculation


## Setup

To run the full 2-D solver, we need:
1) FL distributions
2) strand model
3) count imputation models (RNA and gDNA)
4) a global gDNA density model (the hyperprior)

At setup time:
- finalize our RNA FL and gDNA FL distributions
- fit strand model (BB overdispersion)
- count imputations models remain undefined
- global gdna density hyperprior undefined

### Count imputation model setup

We do not have enough information to fit count imputation models at initialization time. These remain undefined at init

### Global gDNA hyperprior setup

We model the global gDNA density as a single scalar density. However, under hybrid capture, there is no true single global gDNA density. the gDNA density depends on enrichment (targeted by capture probes) or depletion (off-target).

To account for hybrid capture (enrichment / depletion) we can model global gDNA using a single global density scalar and per-node enrichment factors.

- Global gDNA density scalar: the total global gDNA mass divided by the total genome size (bp)
- Enrichment factors: the ratio of local / global gdna density
- Node enrichment factor > 1.0: enriched
- Node enrichment factor < 1.0: depleted

At initialization we have no information. 
- We initialize enrichment factors to 1.0 (uniform).
- We initialize the global gdna density scalar to a numerical floor? **OPEN QUESTION** is this right?

Our initial assumption is uniform gDNA. We set the global gdna density scalar to either ZERO (assume zero gdna) or a numerical floor if there are stability concerns (to be determined)


## Node Initialization

The initialization pass iterates over every node.

***Any node that can be solved directly without imputation is solved***


### Group 1 - Nodes that require no solver (0-D seeds)

1a) Intergenic regions
1b) Intergenic-exon boundaries: { f+=0, f-= 0, fg=1 }

- Fractions: { f+=0, f-=0, fg=1 }
- Uncertainty: { unc+=0, unc-=0, unc_gdna=0 } (propagation sink, no uncertainty)


### Group 2 - Regions that can be solved without imputation (1-D nodes)

2a) Intronic regions (single strand)
2b) Exon <-> Intron boundaries (single strand)
2c) Exon regions (single strand)

Properties: 
- f- = 0 if transcript is on positive strand
- f+ = 0 if transcript is on negative strand

Single-stranded nodes have only 2 active fractions (one RNA and one gDNA). Since we do not have an imputation model at initialization, we solve with the strand model.

To initialize: Apply strand deconvolution ALONE (NO COUNT MODEL, NO IMPUTATION, NO COMBINE). 

Strand deconvolution ALONE will yield fractions and uncertainties based on the strand model (BB overdispersion fit).

For example, a '+' stranded node with U total counts deconvolves `U = Upos + Ugdna`. We compute fractions f+ and fg (f- = 0)

The uncertainty of the deconvolution depends on the strand-specificity parameter and the BB overdispersion.

When ss=0.5 (unstranded), there is no information (infinite uncertainty) and the solution must defer to the weak global hyperprior.

When as `abs(2*(ss - 0.5))` approaches 1.0, we have increasingly reliable deconvolution information. The uncertainty depends on the total counts and the BB overdispersion fit.

- Fractions: { f+, f-=0, fg } where f+ = 0 OR f- = 0
- Uncertainty: { unc+, unc-, unc_gdna } where unc+ OR unc- = 0


### Group 3 - Nodes that cannot be solved by strand deconvolution alone (2-D nodes)

Definition: *any* node with both strands activated: `(ex+ or in+) AND (ex- OR in-)`

These are (strand=AMBIG) and cannot be solved directly by strand deconvolution. They require external information. External sources of information come from:
- Imputation from adjacent nodes via propagation (count imputation model)
- Global gDNA density hyperprior

*At initialization time, 2-D nodes must be initialized to a weak prior*

The default is gDNA with *extremely high uncertainty*
- Fraction: { f+ = 0, f- = 0, fg = 1 }
- Uncertainty: numerically stable MAX uncertainty value



## Model initialization

After node initialization, we have initial fractions/uncertainities for all nodes. We can now build our count imputation models and redefine global gDNA density.

### global gDNA hyperprior

1) Use the initialized results (despite the high uncertainty) to compute the global gDNA density:
- (sum of gdna mass at all nodes) / (sum of node lengths)

2) enrichment ratios

Compute the enrichment ratio for each node:
`(local gdna density) / (global gdna density)`

3) global gdna density hyperprior uncertainty?

**OPEN QUESTION**
- How do we set the uncertainty of the global hyperprior? 

For each node, estimate gDNA mass from the global gDNA density alone: `est_gdna_global = node_gdna_mass * global_gdna_density_mean / node_enrichment_ratio`

Isn't this just going to equal '1.0' everywhere? What am I missing here? We should be able to compute global hyperprior variance/uncertainty somehow.

### Count imputation models

#### gDNA imputation model (var ~ mean fit)
#### RNA imputation model (var ~ mean fit)

The current implementations may be viable, but suspect that they can be improved. 

OPEN ISSUE: we need to audit, understand, and further develop the existing count imputation model fit for both gDNA and RNA.

We have different concepts that need reconcilitation:
- CURRENT: how well do adjacent boundary pairs "agree" on their shared region? This is more or less the current model, I think.

NEW IDEA (TO BE EXPLORED):

Take pairs of adjacent nodes (boundary <-> region). One boundary. One region. Designate one node the 'source' and one node the 'destination', depending on which way signal propagation is going.

Now we would use the current source node to impute the destination node. That means taking the source node's current fraction vector and projecting it onto the destination node. This has passed as a message for the belief propagation system that we've sort of wired already. The destination node has an existing set of fractions in place with uncertainties. So the destination node receives this message from the source node. The message is fractions and uncertainties. The receiving node or the destination receives this message and integrates it with its own fractions and uncertainties.

Forgive me but I am just going to describe this the best I can in plain language.

So the question is, how do we set the uncertainties when we pass a message from source to destination node. Well, the source node already has its own uncertainties, so that needs to be a key aspect. But when the message is passed across, we need to understand what the reliability is or what the kind of communication error is. So the idea of an imputation model is to try to fit the source nodes' imputation estimate with the destination nodes' current fractions. And we can do this with every pair of adjacent nodes. in both directions, really. so we could flip the source and destination nodes and send the message the other way, essentially. So if we have a scatter plot, for example, on the x axis would be the imputed value, and the y axis would be the destination node's actual value.

For each imputation, we have the predicted or estimated fractions and the destination node's actual current fractions. and we would model that we would somehow try to fit the two I think instead of fractions, the scatter plot would sort of be a log log plot with counts on both axes and... or it could be normalized counts with count densities on both axes, but it needs to be log log. And then we would try to fit to that log log plot, and that would sort of give us our message passing or imputation fit. At first I expect the fit to be extremely poor. But, eventually, over iterations, ideally, it would improve.

This is one of the areas where I'm the most uncertain about formalizing the framework. The reason is that for a given pair of nodes, one boundary and one region, I'm not sure exactly whether, like, imputation messages are even needed. For example, an intergenic node is already certain in the adjacent intergenic exon boundary is already certain. So message passing is going to just yield no update, essentially. There's not going to be any belief propagation there. But then if you think of an intergenic to exon boundary node, and then the next node, which is the exon itself, that is valuable to impute. We want to see how well the boundary can predict the exon. We'll have initialized the exon node potentially to an initial genomic DNA value if we have a strand model. If we don't have a strand model, then The exon node will be initialized to a global hyperprior, which would be just genomic DNA with a fraction of one and a very, very high uncertainty. So the imputation fit would be from the boundary onto the exon, and we would use the boundary to impute onto the exon. And on this scatter plot that I'm imagining, they predicted would be on the x axis, and the actual would be on the y axis. We would compute an error between the two. And after accumulating all of these pairs of observations, we could try to fit an imputation model. That model can be used to tell us what the variance or uncertainty is for these imputations. So I believe that at first the imputations will be considered inaccurate because there'll be many in nodes that are initialized to default values. They won't have an accurate genomic DNA estimate. But over iterations, if the genomic DNA estimate improves, the imputation fit will improve. The variance will go down in the uncertainty of belief propagation will go down. So, hopefully, that would lead to convergence.

Now we do need these amputation models for both RNA and genomic DNA. The RNA side comes from spliced fragments. And so a subset of the pairs of nodes will be one boundary that has spliced fragment mass and its associated region. So when there's spliced fragments at a boundary, we would want to impute the RNA mass on the same strand as the spliced fragment mass. And this is an independent and extremely important valuable imputation. So a single imputation between one boundary and one region communicates a genomic DNA estimate from the source node to the destination node and an RNA estimate from the source node to the destination node. and the spliced fragment case is a special case of imputation that requires its own RNA level imputation model. This can be trained in the same way on pairs of nodes that have splice junctions essentially and spliced fragments. But when we impute between nodes, we need to impute the RNA mass and the DNA mass this is essentially the three term fraction communication.


## Iterative full 2D solver

We have now completed initialization. We have everything we need to run the full 2-D solver:
- All nodes have initial fractions and uncertainties.
- We have an initial var ~ mean gDNA imputation model
- We have an initial var ~ mean RNA imputation model
- We have a global gDNA hyperprior model

To solve nodes, we run our belief propagation algorithm (bidirectional sweep).

Each node must update its initial solution with:
- count imputation (gDNA and RNA) from adjacent nodes
- gDNA hyperprior

This is the iterative solver. After each iteration, nodes will have updated fractions/uncertainties.

After each iteration, we re-compute our models. We must update the global gDNA hyperprior (scalar density, enrichment ratios). We refit the count imputation models.


