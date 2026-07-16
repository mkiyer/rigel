Publish a document to 'docs/calibration' describing this new strategy.

Let me explain how calibration SHOULD work. Then, I need you to audit and ensure that calibration is working this way.

At initialization, we have no value for nodes. We initialize them to 100% gDNA at 0 precision by default. These are zero precision nodes. It doesn't matter what the density value is! The precision is zero, so they will move.

Next, we need a message-free self-solve step. This is mandatory. Nodes need to establish themselves and set their internal precision.

# Initialization

## Structure init
- Intergenic regions: pure gDNA, no RNA, infinite precision. Intergenic nodes absorb RNA and communicate gDNA. They are RNA sinks. They are gDNA sources. Unmovable. No messages can pass through because precision is infinite (all messages die at intergenic).
- Intergenic-exon boundaries: same as intergenic regions. These are either transcript start site (TSS) or end sites (TES)
- Single-strand nodes: any node where one of two RNA strands is off. The unexpressed strand gets set to 0 density at infinite precision. Whichever strand is OFF becomes a SINK for RNA on that strand.

## Strand init
- nodes must "self solve" using the strand model, giving them an initial 'tilt' (strand balance) with a defined precision.
- for unstranded data (ss=0.5), precision should be extremely low (~0)
- for strand-specific data (ss >> 0.5), precision should be extremely high

## Boundary init (**NEW**)
- boundaries are special because they directly observe spliced fragments (mature RNA) and unspliced fragments (nascent RNA + gDNA). Nascent RNA is extremely sparse in real biology for unspliced boundary crossing fragments ~ gDNA is a reasonable starting assumption (data can change this later, this is just an initialization assumption)

- Boundaries should then be able to define an initial solution with precision. Should have values for { RNA+, RNA-, and gDNA }. Boundaries need honest internal precision. This is the precision of a DENSITY (counts / length). We need to incorporate honest models following discrete count principles for the boundary initialization. This should be  decided by count (poisson model.. but should incorporate overdispersion.. poisson variance is too generous for precision) over length. Boundary length is quite short, in fact, equal to the FL distribution of RNA (spliced) and gDNA (unspliced).

- If the FL distribution is WIDE, then the density estimate is quite variable. For example, we could have 1 count from a 100bp fragment or 1 count from a 300bp fragment, a 3X variation in density from one count. 

- WE MUST model the precision honestly and this warrants a derivation. 

- Then, boundaries can initialize self-initialize not only with strand, but also with the observed spliced vs unspliced counts. The spliced fragments are already stranded. The strand is fixed by the splice junction, which has a genomic splice motif that is single-stranded (e.g. GT/AG is the canonical motif). So we KNOW the strand of the spliced fragments and then placed them directly into the solved density for the appropriate strand. 

- For the unspliced fragments, IF we have strand-specific data, we can deconvolute nascent RNA from gDNA. Strand deconvolution has a precision associated with the beta-binomial model. Otherwise, we can simply assign everything to gDNA.

- We still need to assign precision to our density estimate, whether it was derived from strand deconvolution or simply count-based in the case of unstranded data.


# **NEW** gDNA prior init (KDE based)

After node initialization, we need to build an initial KDE. This is a new addition. We will build a "pass-0" KDE without even trying to run an initial solve with message propagation.

The solution mirrors our current KDE solution. The current KDE is built after pass-1 (pass-1 is the first solve with message propagation). However, we will built a de novo KDE without our initial message propagation solve. We must do this because we need an initial estimate for our gDNA background hyperprior. The hyperprior is essential to the solve.

Remember our solve is strand balance likelihood, messages, and gDNA hyperprior

Algorithm steps:
1) collect all single-strand exon-intron boundaries
2) extract spliced density (counts / RNA FL dist) and unspliced density (counts / gDNA FL dist)
3) project fraction of unspliced / total density onto adjacent exon region
4) so if single-strand exon-intron boundary has spliced density =4, unspliced density =2, then fraction is 

build a KDE of unspliced density (x axis)




The general node solver is now fully specified and its two theory gates are cleared (count-conservation self-defense + uniqueness). The remaining items before touching production are all sweep-level, not solver-level: change B (count-conservation message integration), D (locked-seam emission), E (the relay fix), and C (prior-free first pass) — plus the still-open #3 transfer-variance model.

The natural next move is to start wiring the §6 solver contract into production (node_geometry statics/init + simplex_logodds ψ) — the solver itself is proven and init-scoped, so it's the lowest-risk first implementation. Or, if you'd rather keep proving theory before any code, I can take on the #3 transfer-variance / enriched-vs-depleted model next (the last undesigned piece). Which do you want?



=======================



Now remember, we use the NPMLE to obtain genomic DNA prior estimates. This is the density value. We also need a method to estimate the density variance. The reason we need this is because message propagation relies on precision, and we don't yet have a strategy to estimate message precision. We have good methods to estimate the precision of observed counts and densities, but when we communicate between adjacent nodes, we also need some sort of message propagation variance. Right now, we model that as the variance of the disagreement between adjacent nodes, but this depends on the density of the observed nodes.

Examples:
- a depleted node pair has high agreement and sits near the depleted gdna floor. variance should be LOW. precision high.
- an enriched node pair has moderate agreement and sits somewhere in the enriched gdna peak. variance should be based on the population of enriched nodes (excluding depleted nodes). variance should be moderate. precision moderate
- a node pair has high disagreement - one node is enriched and one is depleted. the deplete node project to the depleted floor. the high density node projects to an enriched mode. the modes differ. variance HIGH. precision should be very low.


So I think node pair disagreement benefits from a gDNA prior projection function.

We have two nodes A and B

gdna prior (A) = NPMLE_project(A)
gdna prior (B) = NPMLE_project(B)

Do you have ideas on how to do this properly?

- Message precision should be HONEST
- Message precision should not be CONSTANT
- Variance should be derived. It should be a function that depends on density or projection of density onto NPMLE, or something like that?


===========

LEt me articulate the prior setup:

# The prior depends upon "self solved" nodes

"self-solve" means nodes can solve for their own RNA vs gDNA content without a prior and without message propagation.

- With strand-specific data, single-stranded nodes can solve

- Splice junction boundary nodes can ALWAYS solve themselves. Spliced fragments are RNA (strand known). Unspliced fragments are nascent RNA + gDNA. Nascent RNA can be peeled by strand (if available) or density discrepancy (compare to intergenic + intronic regions). 

- "Structural" nodes are frozen fixed and "solved" for eternity.

# What other "solving" can we do before constructing the NPMLE prior?

- if available, use strand-specific  data to solve single-strand nodes (prior free)

- ALWAYS solve single-strand boundary nodes (prior free), just needs background distribution (intergenic + intron regions) to peel nascent RNA from unspliced crossing fragments.

- Limited form of message propagation can be done prior-free to impute exon regions from the adjacent solved boundaries

# What can't be solved without a prior?

- overlapping transcripts in opposite strands create AMBIG nodes. These need message propagation and an accurate prior


# What initialization nodes can we trust?

- Intergenic: pure gDNA, frozen, fixed, infinite precision. Likely fully depleted in hybrid capture.

- Intergenic <-> exon boundaries: pure gDNA, frozen, fixed, infinite precision. Some could be partially enriched in hybrid capture.

- Intronic: gDNA + nascent, SAFE to initialize as 100% gDNA. Strand can peel nascent RNA. Large density discrepancy compared to background gDNA can also peel nascent RNA. Unstranded, we will knowingly and correctly siphon some nascent RNA -> gDNA. Likely depleted by hybrid capture

- Single-stranded Splice junctions (Intron <-> Exon boundaries). Critically important. Gives us COMPOSITION at fixed point. Can project composition onto adjacent exon region with pure arithmetic to initialize -- this is a GOOD first approximation/initialization state and rescue some exon regions. The simple arithmetic composition transfer gives us an initial deconvolution of the adjacent exons, resolving RNA (from the splice junction strand) and DNA. This is the way to allow exon region densities into pass-0 NPMLE.

The unspliced exon <-> intron boundary crossing fragments *can* be nascent RNA + gDNA. Same nascent RNA peel opportunities exist as for introns (strand and density discrepancy).

- Single-stranded Exons with adjacent splice junction boundary nodes. Exons with adjacent splice junction boundaries can make a crude, non-statistical, arithmetic projection of its boundary splice junctions. This is a form of non-statistical message passing.. it assumes the exon has zero precision and completely trusts its adjacent boundaries (single gauss seidel step). An exon may be bounded by ONE OR TWO splice junctions. If two splice junctions, unclear which to trust. Composition can be the average/median/geomean, etc. Unclear. The imputation assumes that the exon has zero precision (trusts nodes completely). 

This IS a prior-free belief propagation step (gauss seidel). Not the full solve. If we do this, we SHOULD give the exon a precision. How? Well, we should have already derived boundary precision. So the exon should inherit the precision of its boundaries. 

It should be implemented as a single-node solve with incoming messages from left boundary and/or right boundary.


=================


Let's review what I see are the FOUR missing pieces (there may be more, this is my mental map, correct me if I'm wrong)

# 1) No RNA prior

It was absolutely correct to remove the "RNA message on log(1−f_g)". We are on track.

We need to establish the CORRECT derivation:
`ψ = strand + logP_g + logP_r`

This is the answer.

So we need `logP_r`.

That is a missing piece.

# 2) Respect message interdependence

We are solving for a PIE CHART with 3 pieces: RNA+, RNA-, and gDNA. Changing one value FORCES the others to change.

We MUST have our solver respect this.

What happens if RNA precision = 0, but DNA precision = HIGH.

Then we receive a weak message with nonzero RNA precision?

It will MOVE the DNA value!! This is wrong. The DNA precision is HIGH.. the solver should not be allowed to change it when a weak RNA message arrives.

This is the 1-D version (1-DOF)

AMBIG nodes have three pieces of the pie to solve for. There are 2-DOF. We must respect this.

# 3) Spliced mass needs to participate!

When we pass messages boundary -> region, what happens to the spliced mass?

The spliced mass combines with unspliced mass into a single message. The message is received by the region node and used by the solver.

Example of a splice junction on the POS (+) strand:

Boundary:
- spliced density (by definition this is RNA+) = 10
- unspliced densities { RNA+ = 1 (nascent RNA), RNA- = 0, gDNA = 10 }

How does the message get passed?

the total RNA density needs to be aggregated into the message. There may be a fragment length correction that i'm skipping here, but the concept is that the splice + unspliced mass = total rna mass

So the message will have { RNA+ = 10 (spliced) + 1 (unspliced), RNA- = 0, gDNA = 10 }, along with their precisions.

The REGION will receive this message. It receives the total RNA density from the sender. That should solve the wiring issue

# 4) Message precision is still broken

This still needs to get fixed. But it's hard to wire in proper message precision until we fix everything else.

What should be the path forward here?

I think we need to STRIP the wrong variance term out of messages while we develop this set of fixes.

If we MUTE our messages, then we won't see the full cascade effect of the message propagation.

So for the time being, while we have simulated synthetic data with ground truth, I would remove the shackles from the message precision. Allow messages to run wild (let their precisions explode if need be). 

Unconstrained message precision will allow us to safely develop the other fixes in our calibration system. If we have bugs or issues, we will see them propagate in messages.

I think muting / constraining our message should be one of the last steps in development.






