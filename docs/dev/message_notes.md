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



















































