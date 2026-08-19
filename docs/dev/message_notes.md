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














































