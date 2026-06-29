
# Goal of calibration

The goal of calibration is to deconvolute unspliced fragments into gDNA and RNA.

Given that RNA is single-stranded, we must deconvolute unspliced fragments into RNA+, RNA-, and gDNA.

The calibration result forms the bayesian prior for the EM algorithm. We need to tell the EM how much gDNA to expect. We have proven that accurately estimating the prior leads to accurate abundance estimation.



# Training nodes

Some nodes can be solved without additional information:
- Seed nodes (fixed by transcript annotation structure)
- Single-strand nodes (with strand-specific data)
- **Splice junction boundaries** can be partially solved because the spliced fragments at the boundary are pure RNA. The unspliced fragments must be gDNA (or nascent RNA, assumed to be sparse/low in most protocols).

Seed nodes:
- Intergenic regions -> pure gDNA (assumed)
- Intergenic <-> exon boundaries -> pure gDNA (assumed)
- Intron regions -> gDNA (and sparse nascent RNA)
- Intron <-> exon boundaries (gDNA and sparse nascent RNA)

## Exons

What about exons?

In a single sample, many transcripts are expressed, and many transcripts are not expressed. There are ~500,000 annotated transcripts in the human genome. The majority of transcripts are not expressed, or have very low abundance. A small minority of transcripts consume the vast majority of the RNA-seq counts.

If a transcript is not expressed, its exons will have zero RNA. All unspliced fragments can be assigned to gDNA.

If a transcript is expressed, its exons will have nonzero RNA. The unspliced fragments are a mixture of gDNA and RNA.

So:
- Unexpressed exon -> gDNA
- Expressed exon -> gDNA + RNA

The problem is, we don't *know* which exons are unexpressed. We have to solve for this.

Single-strand exons deconvolve directly into RNA + gDNA. Strand-specific data solves it.



# Overview of hybrid capture

When we have hybrid capture data, we have probes that target  nucleic acid from specific genomic regions.

This causes enrichment of those regions, and the depletion of the off-target regions.

In the human genome, there are now ~75,000 genes and >500,000 transcripts.

Hybrid capture tends to target at most 20,000 genes. The majority of the annotated transcripts will be off-target.

We must assume that hybrid capture thus targets a subset of the annotated transcripts.

The targeted subset will be enriched roughly 10-1000X from their baseline values. The off target subset will be relatively depleted.


# Non-capture versus capture

With hybrid capture, we have a subset of enriched nodes (on target). The remaining nodes are relatively depleted.

Conceptually it helps to divide nodes into one of four categories based on RNA abundance (zero, non-zero) and hybrid capture (on-target, off-target).

Here we wish to model hybrid capture by an enrichment ratio. This is because hybrid capture tends to be on a spectrum. Some probes capture more strongly than others. Not all regions capture the same. 


# Detecting enrichment by monitoring gDNA

RNA-seq already has enormous dynamic range of >10,000 fold. Therefore it is challenging to determine captured vs off-target RNA from a single sample. Hybrid capture *could* amplify the already huge dynamic range even further, but it is more likely to shift the relative abundances of transcripts than it is to create bimodal (enriched vs off-target) behavior.

This makes the question, "do we have hybrid capture data?" difficult to answer from RNA alone.

Genomic DNA is another story. With an approximately uniform genomic distribution (plus local variation by GC content and other sequencing factors), it becomes plausible to detect enriched regions by monitoring gDNA levels.


# Belief propagation / message passing

We rely on imputation using belief propagation to solve any nodes that cannot be directly solved.


# Method enhancements and setup


## Global gDNA minimum (not mean)

Currently, our algorithm lacks a global gDNA density minimum. We measure and use a global gDNA density *mean*. This works with non-capture data. With capture data, we have a subset of highly enriched nodes and a subset of depleted nodes. The mean averages the two.

1) intergenic gDNA density mean - this forms our true gDNA density minimum. This is the sum of unspliced fragments over the total length of intergenic regions. Intergenic regions are assumed to be relatively depleted under hybrid capture.

2) intronic gDNA mean - introns should closely approximate the intergenic gDNA density mean. Nascent RNA can be present in introns but typically at low levels, and quite sparse genome-wide. In real data, some transcripts have measurable levels of nascent RNA, but these tend to be quite rare.

Together, the intergenic and intronic nodes can be used to estimate the gDNA density minimum in the dataset (the depleted, off-target level of gDNA). The amount of gDNA should not go below this.

This is the "floor" that model fitting must shrink to.


## Attempt to "Self-solve" nodes at initialization time

During init/setup, all unsolved nodes should solve themselves using 1) structural signatures, 2) the strand model.

Structural signatures (such as intergenic) translate to infinitely precise densities.

The strand model may or may not be usable.  With unstranded data, this should yield a GARBAGE solution with ZERO precision. With excellent strand-specific data, this can yield a highly precise solution.


## Build belief propagation models for imputation

Consider an unsolved internal exon region R and its intron <-> exon boundaries B1 and B2.

Can B1 or B2 predict R?

Let's think through the scenarios.

### Non-capture

Unspliced crossing density of B1 is gDNA (+ low nascent RNA)
Spliced crossing density of B1 is RNA

There is no hybrid capture, and so the density of R should be approximately equal to the densities of B1 and B2.

We should therefore be able to DIRECTLY IMPUTE the gDNA and RNA densities of B1 onto R.

Our current modeling fits log(R) ~ log(B1 or B2)

The residuals form the variance. 

**The error tells us how to set the precision** -- this is something we are still working

Our latest results suggest that this is working well!

### Capture

Work in progress.

