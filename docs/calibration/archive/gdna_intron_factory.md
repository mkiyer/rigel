
# the problem

We can trust the intergenic regions are pure DNA, and we build an intergenic gDNA distribution using these regions. This is our BACKGROUND gDNA distribution

Nascent RNA refers to active transcription that has not been fully processed. Transcription produces RNA that spans an entire transcript from start to end. Splicing machinery cleave the introns, a poly-A tail is added, and a 5' cap is added. This is the transcriptional processing that happens. Nascent RNA refers to the unprocessed immature transcript.

Detection of nascent RNA is possible by measuring intronic fragments. In real rna-seq data, nascent RNA is very sparse and usually not substantial. However, in rare transcripts or certain library prep procedures, nascent RNA can be readily detected.

Rigel models nascent RNA, but the model is not yet complete. We wish to extend our modeling of introns, where nascent RNA can be unambiguously separated from mature RNA (Exons have mature RNA, nascent RNA, and genomic DNA).

# the idea

Introns harbor genomic DNA at background levels, PLUS nascent RNA. We can therefore detect gDNA within introns by comparing the unspliced fragment abundance to the intergenic fragment abundance.

The goal is to deconvolve introns into genomic DNA and nascent RNA components.

- in real data, there is scarce nascent RNA
- when nascent RNA is present it must be above background gDNA levels
- we must detect gDNA within introns. this is an especially critical need for unstranded data

# concept plan

1) model intergenic gDNA (the background)

- model background gDNA levels using intergenic regions

2) model intronic gDNA (gDNA + nascent RNA)

- model intronic regions

3) deconvolve individual introns

the unspliced fragments within an intron should be modeled as being drawn from a mixture of two distributions: the gdna background distribution (intergenic regions) and the intronic distribution (gdna + nascent RNA).

genomic DNA lives EVERYWHERE. Introns therefore harbor both genomic DNA AND nascent RNA.

So the goal of the method is, given an INTRON region, deconvolve its total counts into nascent RNA + gDNA.

We can do this based on likelihood. 

We may be able to use our existing density NPMLE for this:
- intergenic regions -> density NPMLE
- intronic regions -> density NPMLE

then for each intron, ask probability of being sampled from intergenic vs intronic?

But remember that intronic total = intergenic + nascent RNA.

So given an observed density, how much is gDNA? At what confidence level?


# goal

1) understand the problem
2) derive candidate solution ideas
3) model intergenic and intronic distributions
4) prototype (deconvolution of introns)


# constraints / needs

- deconvolution of introns into gDNA + nascent RNA must have a precision associated with it. the precision must be honest, honoring the discrete count over a discrete length nature of the data

- the deconvolution should produce a CONFIDENT (high) precision for gDNA in most cases

- I would expect it to be imprecise about peeling off nascent RNA, as the intergenic vs intronic distributions are mostly overlapping


# open questions / problems

1) do we need to model both intergenic and intronic distributions to perform the deconvolution? the global pattern of intronic abundance does not inform an individual region. individual introns may harbor nascent RNA and that does not really correlate with a global pattern. If we decide to peel off nascent RNA based on the intergenic distribution alone, then we need some kind of magic number "cutoff"? how do we assign precision to the intronic fragments after deconvolution?

2) assigning strand!

RNA is single-stranded. If we peel of nascent RNA from an intron, we don't know which strand to assign it to. It's better to think of it as peeling off gDNA! 

"We have a pool of unspliced fragments, how many fragments are confidently gDNA?"

We can do this! but we can't assign strand to the remaining nascent RNA. If it's a single-stranded intron, then the strand can be inferred. If it's an ambiguous intron, then the undesignated nascent RNA must be resolve by the full solver.

Similar to the idea of strand "tilt" in an ambiguous region, intronic gDNA deconvolution informs us of the "tilt" in an ambiguous (both strands active) region. It can tell us "this intron is likely 95% gDNA" with an associated precision, but it can't tell us what happens to the remaining 5%.

In a single-stranded intron, "this intron is 95% gDNA" implies "the rest is RNA on the inferred strand"


# Summary of this "nascent factory" idea

- use intergenic vs intronic abundance to "peel" confident gDNA within introns
- this should manifest as a gDNA DENSITY with a gDNA PRECISION -- in our new precision model we have 1-DOF (for single-stranded nodes) and 2-DOF (for both-stranded ambiguous nodes). the belief must reflect the mode and precision of intronic gDNA


# Extending to all nodes

## idea

Using the intergenic gdna background level can apply to all nodes, and not merely introns.

Your review highlighted that the approach can be extended to exons (we can peel gDNA out of exons in the same fashion). 

**The trouble is hybrid capture**

Hybrid capture enriches all nucleic acid (both RNA and gDNA) above background levels. The intergenic distribution only measures background gDNA levels. So while it is not wrong to peel confident gDNA within exons using the exact same approach, much of the total density within exons will be unaccounted for, because all we can confidently peel is the depleted floor density. 

In otherwords the intergenic gdna factory implementation can be applied globally to all non-intergenic nodes, but will underestimate gdna under hybrid capture conditions.

If we apply our gdna factory to exons, we must be sure to derive the correct precision term.

For example:
- suppose intergenic background gdna has mean, variance
- exons and intron-exon boundaries are enriched by a factor of 1000X over the background level (background is depleted)

If we apply our gdna factory to enriched exons and intron-exon boundaries, we can successfully peel gDNA.

But what happens to precision?

- A exon or exon-intron boundary node might have total density=30
- Our intergenic background may be mean=0.01 (with variance)

- We can peel off the intergenic background (approx 0.01), leaving majority of density (30 - 0.01) unaccounted for. The method can say, "i'm confident that 0.01 is gDNA", but it cannot say anything about the residual 29.99 density. The residual may still be gDNA or RNA.

So when there is a HUGE discrepancy between total densities, the gdna factory must remain HONEST about its capabilities. It can still work, but cannot result in a confident call about the total density composition of the enriched nodes.

How does the method work now? Can it be extended to ALL nodes safely, without making the enriched nodes overconfident about their final solution?



