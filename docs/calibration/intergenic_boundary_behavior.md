# what is STORED in a node's current belief?

a node's current belief can be stored as a fraction (composition) summing to 1.0, OR as raw counts. The two are interconvertable easily.

we do NOT store density. density is computed from counts and the appropriate effective length. RNA and gDNA have different effective lengths, and so as the relative composition changes, the total density changes. 

so the current belief should be stored as fractions OR counts, along with two precision values. The first is the precision of the gDNA level. The second is the precision of the "tilt" or strand balance between RNA+ and RNA-. ***NOTE: verify this is true**

the 'gate' (gdna always on, rna+ on/off, rna- on/off) determines whether the solver needs 1-DOF or 2-DOF. When RNA+ and RNA- are both gated on, the 2-DOF solver is required.


# principles of message passing

## total density is not conserved

- a node's counts are fixed
- gDNA and RNA may have different fragment length (FL) distributions
- therefore, as the composition of a node changes, the total density can change


## what should be passed in a message?

- my intuition says that messages should be *passed* as density values. the densities are real measurements of absolute levels. this is a common currency that travels between nodes.

- we pass 3 densities: { RNA+, RNA-, gDNA }

- we also have to pass some form of CONFIDENCE measure (the precision) about each of the densities. we might be quite confident about some of the densities but know nothing about the others.

- so we pass 3 precisions { pr_rna+, pr_rna-, pr_gdna}. is that right?

- then, we need to tell the recipients which components we have the ability to measure. for example, intergenic regions can only measure gDNA. Single-stranded nodes can only measure RNA on ONE strand. The message content must be the source's signature: { rna+ on/off, rna- on/off, gdna on/off }

## full message form:

A message contains 3 tuples with (density, precision, gate):

- RNA+ (density, precision, gate)
- RNA- (density, precision, gate)
- gDNA (density, precision, gate=TRUE)

** DO WE AGREE ON THIS FORMAT **

- if we agree, this solidifies how messages should be constructed and SENT. 

- The gate is constant, and always True for gDNA, It doesn't need to be stored in the message.


## implementation

- we must agree on a message format
- message construction and emission should be the first task.

The *implementation* may differ from the principles outlined above. In particular, the "gate" term is a structural constant and never changes. It is necessary to know whether a term is on/off, but this doesn't necessarily need to be encoded in the message. It can be read from a constant vector for efficiency purposes.


## should we gate message emission?

- my intuition is 'NO'. the source node does not gate message emission. 

- the source node always emit messages

- however, we may later chose to add gating logic *ONLY for performance reasons* if the results are bit-identical. for example, messages with ZERO precision across all components need not be sent, because they would be ignored by the recipient.


## how do we interpret messages?

a destination (recipient) node receives a message. the destination must then interpret the densities in its own frame.

## principle of non-uniformity

we ASSUME that two adjacent nodes are NON-UNIFORM in their density distribution.

**source density != destination density**
**DENSITY IS NOT CONSERVED ACROSS NODES**

## how do we handle non-uniformity?

when hybrid creates enrichment/depletion "cliffs", assume that the composition of nucleic acid does not change across the cliff, only the total density

this assumption is not always true, but it is the only assumption we can make given our constraints


## principle of composition

density composition is conserved between adjacent nodes with *equal gates*

the 'gate' tells whether a component is active (on/off).

if source nodes 3-term gate vector == dest node 3-term gate vector, then density composition is conserved.

if source nodes 3-term gate vector != dest node 3-term gate vector, then density composition is not necessarily conserved.

## computing the composition log-odds shift

To transfer composition (ONLY when the 3-term gate vector is EQUAL between two nodes)

we use the composition shift — imputes the dst's fractions from the source composition, with an eff-length frame correction so the capture enrichment cancels (cliff-invariant).

NOTE: this arithmetic is complex and the derivation needs to be verified.

Derivation (needs verification):

src_f_g = src_rho_g / src_rho_tot
src_f_rp = src_rho_rp / src_rho_tot
src_f_rn = src_rho_rn / src_rho_tot

dst_rho_g = src_f_g * dst_count_g / dst_eff_len
dst_rho_rp = src_f_rp * dst_count_rp / dst_eff_len
dst_rho_rn = src_f_rn * dst_count_rn / dst_eff_len

This transformation is done BY THE RECIPIENT!


## how do we receive messages when the gate is not equivalent?

simplest example is intergenic-exon boundary (pure DNA) -> exon (DNA + RNA)

Sent (from intergenic-exon boundary):
- gdna = (rho, pr, gate=True)
- rp = (rho=0, pr=infinity, False)
- rn = (rho=0, pr=infinity, False)

Received (exon region):

Say the exon is on the + strand for this example.
We have *unequal active components / gates*. So we cannot directly use composition shift!

How do we utilize this information? The only active term is gDNA (gDNA is always gated on). 

The recipient now knows that there is measurable gDNA. The fractional composition cannot be computed, because the RNA components are not known.

Do we let the message impute composition?

**This is our dilemma**

We cannot impute composition when the activity gate is not equal. How should we handle this? We need an elegant solution here.

Example boundary sends density = 2
Exon receives density = 2 but has total density = 200
There is 100X discrepancy between boundary and exon.

The difference can be explained by either RNA or gDNA
- exon region could have gdna=2, rna=198
- or exon region could have gdna=200, rna=0

We have no way to tell the difference without an initial gdna prior.

- So what needs to happen?

KEY FACT:
- this region is UNIDENTIFIABLE without a gDNA prior.
- If we update the exon node belief, we will contaminate the solution with noise.
- we must WAIT until the gdna prior is fit before solving.

We maintain two precision values for the node belief.

The first is the precision of the gDNA level. This term includes the boundary counts, the transfer variance, and perhaps other terms. It must be very weak.

The second is the precision of the tilt between rna+ and rna-. The incoming message is gated OFF for rna-, so the tilt can be set to 100% RNA+ with infinite precision.

This leaves a 1-DOF problem with RNA vs gDNA in the exon. From this message ALONE, the solution will be 100% gDNA, at very weak precision.


### single-exon transcripts

when we have unstranded data, single-exon transcripts are fundamentally unidentifiable. without strand information and without splice junction information, there is no direct measurement of RNA signal.

this problem can be partially extrapolated to transcript start sites (TSS) and transcript end sites (TES) even when we have multi-exonic genes.

how *should* we handle TSS and TES for single-exon transcripts?

If we handle it as described above, then single-exon transcripts may only receive DNA messages from their boundaries. they will be set to 100% gDNA with a weak precision (determined only by their boundaries).

I would argue that this is CORRECT behavior, in the absence of a gDNA prior.

For pass-0 (prior-free), this is the correct default.

Single-exon transcript should be WITHHELD from training a gDNA prior, because of this known limitation. A fundamental identifiability problem.


# how does the gDNA prior work?

Imagine a single-exon transcript that over-called its initial solve to 100% gDNA (weak precision determined by the boundary crossing messages).

Let's say now that we train a gdna hyperprior with two modes, density = 2 (depleted) and density = 20 (enriched)

Let's say the total density of the single exon is 30.

We fit the gDNA hyperprior (single-exon transcripts are held out), and now we re-solve.

On the repeat solve, we have incoming messages again from both boundaries.

left boundary (density 2) -> exon <- right boundary (density 4)

The node was UNIDENTIFIABLE in pass-0 (no prior).

NOW, the fit gDNA prior will contribute to the solution.

src_f_g = 1.0 (intergenic-exon boundary is pure gDNA)
src_f_rp = 0
src_f_rn = 0

dst_rho_g = 1.0 * dst_count_g / dst_eff_len
dst_rho_rp = 0 * dst_count_rp / dst_eff_len
dst_rho_rn = 0 * dst_count_rn / dst_eff_len

The messages say, "100% gDNA" using composition invariance, but weak precision.

What about the prior? Prior has depleted mode at density=2, enriched mode at density=20

This exon region node has density = 30. So what happens if we project density = 30 onto the gDNA prior now? 

The prior to pull us back to a density of 20 (its enriched mode)!

This says, "i know you see density of 30, but gdna can only account for 20". 

Solution should then be gDNA ~20, RNA=10. Correct.


# Summary

- During pass-0 (prior-free) it is dangerous to assign precision > 0 to nodes that are unidentifiable. Any unwarranted precision will cascade, preventing nodes from moving to the correct belief state based on the prior

- The composition shift transformation can overestimate and underestimate DNA across unequal gates. We need to carefully derive the arithmetic for the desired composition shift.



# Boundary <-> Exon cases


## Intergenic-exon -> Exon (covered above)

- intergenic-exon is pure gDNA
- exon is DNA + RNA
- gate inequality, message cannot be received by recipient unless gDNA prior available

## Exon -> Intergenic-Exon

- message dies (intergenic-exon is a pure SINK)


## Intron -> Intron-Exon
## Intron-Exon -> Intron

- The intron self-solves to gDNA + nascent RNA during initialization.
- The intron-exon boundary is also gDNA + nascent RNA
- These are gate-compatible, composition log-odds shift applies

**Messages can be directly passed and interpreted from Intron -> Intron-Exon**

No change to current method

Same paradigm applies bidirectionally Intron <-> Intron-Exon


## Intron-Exon boundary -> exon region

## 1) do we have splice junction (yes/no)

- is boundary a splice junction? if yes, is the splice direction in the same direction as the source -> destination of the message? if yes, then the boundary emits spliced density (it is a splice junction source)

### 2) message content

- unspliced gDNA + nascent RNA (derived from intron)
- spliced RNA (what is the precision?)

### 3) emit message

compute densities:
- spliced_density = spliced mass / splice_junction_eff_len
- unspliced_rnapos_density = f_rna_pos * tot_unspliced_counts / rna_eff_len
- gdna_density = f_gdna * tot_unspliced_counts / gdna_eff_len

modes:
- RNA+ = spliced_density + unspliced_rnapos_density
- RNA- = 0 (for pos strand transcript)
- gDNA = gdna_density

#### precisions:

RNA+: *this is the challenge*

RNA is directly measured at the boundary. It should have "count precision", quite high.

Nascent RNA is deconvolved/imputed. It has lower precision because it's levels are based on likelihood of DNA vs RNA in the intron.

So we have two values with two precisions!

- spliced_rna_density has count precision (directly measured)
- nascent_rna has separate precision (based on strand balance and deconvolution precision)

How do these precisions combine?

If there is 0 spliced mass, then we have unspliced RNA alone
If we have 0 unspliced RNA, then we have spliced alone
But if we have both? Do we average the precisions? 

If we have high spliced density and low nascent RNA, then the precision will approach the precision of the spliced density. If we have low spliced density and high nascent RNA, then the precision approaches the precision of the nascent RNA. 

So somehow the precisions merge. How should this be done?

**this is an open question**

otherwise, this should address the message emission


### 4) receive message

- exon region receives message with RNA (one strand) and gDNA densities and precisions

- theoretically, we should have equal gating between intron-exon boundary and exon region. same active components.

- composition transfer should hold here. the exon should be able to perform composition shift using the densities it receives in the message.


## Exon region -> Intron-exon boundary

### compose message:

- Exon harbors unspliced fragments only
- Message composition is straightforward! Just unspliced densities and precisions.

### receive message:

Boundary receives message

1) is this a splice junction boundary? if yes, is the splice junction on the correct side of the message direction? if yes, this is a splice junction sink (absorb).

2) boundary <-> region gates should be compatible, same components are active

3) boundary must absorb the splice junction

- message rna density must be split across splice junction and unspliced fragments, this is the complexity.

- if this is an enrichment cliff, does the boundary composition change? do we preserve composition? this is an open question.

- the simplest solution, is to have boundary split the RNA message into spliced and unspliced in its own frame.

Two equations:
- exon message RNA density = spliced rna + unspliced rna. 
- boundary total density = spliced rna + unspliced rna + gdna

but boundary total density != exon message density!

So how does the boundary node solve? we cannot simply subtract the spliced rna from the incoming message because of the enrichment cliff. 

The intron-exon boundary integrates two messages:
- message from intron (depleted)
- message from exon (may be enriched)

This is one of the most difficult nodes to solve.

This requires a careful derivation. What do we already have?

- We CAN solve for gDNA and nascent RNA (unspliced) using composition shift directly from the neighboring intron, that respects the enrichment cliff.

- We can't use composition shift from the exon message until we have a gDNA estimate for the node. The gDNA estimate is critical. The best source of the gDNA estimate could be from the exon itself, once the gDNA hyperprior is known. The messages from the intron may be sparse, suggesting that gDNA is scared (depleted mode). It will be an imprecise gDNA estimate at best.

- A gDNA estimate from the exon message will be more accurate. Composition shift can be used for gDNA assuming gates are compatible. 

** IS THAT A NEW THEOREM? THAT COMPOSITION SHIFT CAN ALWAYS BE USED FOR GDNA ** If so, that is a major finding. 

exon message:
- dst_gdna_counts = src_gdna_density * dst_gdna_eff_len



- Absorbing splice density using composition shift requires a gDNA estimate. 


- handle gdna messages (both sides). how does simplex solver 
integrate two messages (forward/backward) currently? average the densities? average the precisions?




- spliced composition is RNA only
- unspliced composition is (RNA + DNA)

rho_spliced = rho_src_rna












## Intergenic: nothing, intergenic is solved by definition
## intergenic-exon boundary: defined, solved
- include pass-0

## single-stranded intron
- ALWAYS solvable, "self solve" using strand AND/OR *NEW* density deconvolute (vs intergenic)
- include pass-0

## ambiguous (both strand) intron
- if strand-specific data is available, can "self-solve" using BOTH density deconvolution (provides gDNA estimate) and strand tilt, making the node independently solvable.
- if unstranded, NOT directly self-solvable (no tilt information), but can be solved by message propagation (provides RNA "tilt") without a gDNA prior, because gDNA is by definition depleted and estimated from the density alone.
- exclude pass-0

## single-stranded intron-exon boundary
- may be enriched by hybrid capture
- if strand-specific data, can "self solve"
- if unstranded data, can solve by message propagation
- include pass-0

## ambiguous intron-exon boundary
- this is a boundary that contains introns in both directions and an exon
- if strand-specific, needs message propagation (gdna estimate) + strand -> solvable pass-0
- unstranded: message propagation + strand (solvable)
- gdna hyperprior can stabilize


## single-stranded exon
- may be enriched by hybrid capture
- solvable pass-0
- if strand-specific, can "self solve"
- if unstranded, needs message propagation


## ambiguous exon
- may be enriched by hybrid capture
- not solvable pass-0
- if strand-specific, needs strand-balance + messages +/- gdna hyperprior
- if unstranded, needs strand-balance + messages + gdna hyperprior

## exon-exon boundary
- depends on single-stranded vs ambig

## single-exon transcript (TSS and TES)
- strand-specific: self solve
- unstranded: requires messages + gdna hyperprior






