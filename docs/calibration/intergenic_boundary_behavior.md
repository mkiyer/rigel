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


## Exon -> Intergenic-Exon

- message dies (intergenic-exon is a pure SINK)




