Can you explain what "mechanism 1" is for the new message propagation policy?

I have an idea for how to approach the new message propagation policy.

I want to handle one node type at a time:
1) exon regions
2) intron|exon boundaries (an intron's flanking boundaries)
3) intron regions

For each node type, we need to define:
1) what message does this node SEND? what is the message precision?
2) how does the node SOLVE with its own belief, two incoming messages, and gdna prior?

AS you know, SENDING (propagate()) is different than SOLVING.

I want to start with exons, because they are the most difficult.

Here are my questions. When our strand model is used to solve an exon, what precision does the exon get? I want to see the raw numbers. Show the an e

=========

2026-09-01 message policy development

# Initialization

At initialization we solve nodes:
- intergenic (structural)
- intronic (density model, compared to intergenic background)
- single-stranded (strand model, only when strand-specific data is available)

# Pass-0 solve

On the test chromosome, what isn't solved?

- intron|exon boundaries
- exon regions

That's all. So we need messages to solve these two nodes.

## intron|exon boundaries

Let's start here.

An intron|exon boundary has:
- one flank is an EXON (exon bit set)
- one flank is an INTRON (NO exon bit set, INTRON bit(s) set)

There are TWO incoming messages:
- EXON message that comes FROM the exon source TO the intron|exon boundary destination
- INTRON message that comes FROM the intron TO the intron|exon boundary

Let's begin with a simple solver:
- ZERO the EXON message (make the exon message do nothing)
- ONLY use the INTRON message (allow the intron message to influence the solve)

How do we use the intron message?
- The intron and the intron|exon boundary SHARE COMPOSITION. By "share composition" we mean "THE EXACT SAME TRANSCRIPTS OVERLAP BOTH THE INTRON AND THE INTRON-EXON BOUNDARY". This is a key definition.

*node pairs that share composition can be solved by composition transfer*

In otherwords, if the intron is 90% gDNA and 10% RNA, we directly transfer this composition to the intron|exon boundary.

So how do we proceed?

Again, to start, we will ZERO the exon (source) -> intron|exon boundary message. This is to simplify the solve. We will only use the INTRON -> intron|exon boundary message.

We solve by direct composition transfer.

The main question becomes -- How do we assign precision? This becomes the key question and is the problem to be addressed first.

Composition transfer needs an HONEST PRECISION mechanism.

So in this turn, we hold everything else constant. We apply composition transfer from intron -> intron|exon boundary. We allow that to solve the intron|exon boundary.

My plan would be to assign a miniscule precision to start just to see how that single message affects the solution. Assign a tiny tiny precision (1e-10), something nonzero.

Again we have the exon -> intron|Exon boundary message OFF for now, so we isolate only one message, one kind of solve (composition transfer), and can derive the precision of that one-hop message.

==================

## exon regions

An exon REGION has two flanking BOUNDARIES. Message propagation critically depends on these boundaries.

Solving with message propagation means the EXON region needs to determine how to handle incoming messages. At SOLVE time, the EXON region has TWO incoming messages.

Both messages come from the flanking BOUNDARIES.

Let's discuss the characteristics of the boundaries:
- transcript start or transcript end flag (TSS or TES): present or absent
- splice junction on the same strand: present or absent

Which boundaries can use COMPOSITION TRANSFER?
- no TSS or TES flag
- same strands active on both sides of boundary


Boundary -> Exon message:
- spliced fragments
- unspliced fragments


An incoming message carries information about SPLICED FRAGMENTS (pure RNA) and unspliced fragments (may be RNA or gDNA).

Spliced fragments enter the EXON via a splice junction. Splice junctions are single-stranded. Incoming spliced fragments SPLICE IN to a boundary. The term "splice-in" means that spliced fragments enter an exon region.

So an incoming message carries:
- spliced fragments (single stranded RNA)
- unspliced fragments (RNA+, RNA-, gDNA)

An EXON region by definition ONLY has unspliced fragments. So the incoming message conveys both spliced RNA and unspliced RNA. The Spliced RNA and Unspliced RNA is aggregated together and imputes the EXON's RNA. The unspliced gDNA imputes the exon region's DNA.


### sending a spliced message

In previous sessions, we established the the SPLICED-IN message is the aggregation of all the incoming splice junctions.

For example:

TA+ exons (1000, 2000), (10000, 11000)
TB+ exons (5000, 6000), (10000, 11000)

Exon(10000, 11000) has TWO incoming splice junctions, PLUS unspliced RNA. The Spliced-in RNA is the aggregate RNA from each one of the incoming splice junctions. Additional unspliced RNA may come from introns OR exons.

Provided that the boundary at position=10000 has no TSS/TES flag AND there is no opposite stranded transition (- strand), we can use composition transfer between the BOUNDARY at position 10000 and the exon region(10000, 11000).

spliced fragments are STATIC and NOT "solved". They are fixed. They add signal to existing messages.


The first goal is to design the architecture for the message propagation INTO exons. The message must incorporate the SPLICED-IN RNA, the UNSPLICED RNA, and the gDNA (unspliced by definition). The message honest precision must be derived.

This is our first task.


# intergenic|exon boundaries

we need to briefly discuss and confirm that intergenic|exon boundaries are certified pure gDNA structurally. They are always 100% gDNA by definition. They are not solved.


## exon regions with terminus boundaries

There are many situations where an incoming message does not have composition-compatible information. For these messages, we cannot employ the composition transfer paradigm.

We need to dissect the individual situations where we have incompatible messages, and address them one by one.

Let's start with what is probably the simplest: the intergenic-exon boundary.

### intergenic|exon boundary <-> exon message

We are solving for EXON regions. When an EXON region is bordered by an intergenic|exon boundary, we cannot use composition transfer.

Characteristics of the Intergenic|exon boundary:
- Only gDNA (unspliced) crosses the boundary
- The boundary MUST also be a terminus (by definition)
- intergenic regions are currently defined as pure gDNA (structurally, due to complete lack of annotation overlap)

Genomic DNA crosses contiguously across the boundary as unspliced fragments.

Hybrid capture probes placed near the intergenic|exon boundary can partially enrich these fragments. Typically, the probe will not fully enrich the boundary crossing fragments due to partial overlap.

So how to we impute from the intergenic|exon boundaries onto exons?

- We can directly impute the gDNA abundance itself from the boundary onto the exon.
- The precision of the imputation depends on the same precision determinants (count, effective length, etc) similar to the other cases.


### how does an exon reconcile TWO messages?

when we add intergenic|exon boundary messages, we will have more exons that have to reconcile two message sources.

If our structure is:
<intergenic> | <exon> | <intron> | <...>

The <exon> receives messages from the intergenic boundary and the intron boundary.

Prior to the initial solve, in the unstranded case, the exon has no 'own' belief. It relies solely on messages.

The intergenic message provides a gDNA level, but the message cannot provide an RNA estimate.

The intronic message provides a composition (gDNA and RNA).

At SOLVE time, the exon has both messages.

There is an option to *BORROW* composition information if it exists in either one of the two messages.

The intron boundary message carries composition information. The intergenic boundary message carries only gDNA abundance. 

**How can these be reconciled?**

This is the crux of the question.

Solving this is the goal of this next sequence of turns.


### single-stranded exon|exon boundary

An exon|exon boundary implies that there are EXONIC regions on both sides. There can also be an intron on one or both sides too.

Example:
TA+ (1000, 2000), (20000, 21000)
TB+ (1000, 2000), (10000, 11000), (20000, 21000)
TC+ (500, 22000)

The simplest way to create lots of exon|exon this is to cover everything with one long transcript span (TC+). Now every boundary is exon|exon. 1000, 2000, 10000, 11000, and 20000, and 21000 become exon|exon boundaries.

We can consider each case and make sure we have the right model. Let's do this one case at a time, and build the model one case at a time.

Let's start with the Boundary at 1000:
- left region is exon (5000,1000)
- right region is exon(1000,2000)
- terminus present (TSS+)
- if message propagation is FORWARD (left -> right), we ADD a new transcript (transcript gain)
- if message propagation is REVERSE (right -> left), we LOSE a transcript (transcript loss)

We have TWO problems to address: (1) message propagation, and (2) solving.

Let's solve for boundary 1000:
- two incoming messages, left is (500,1000), right is (1000,2000)
- boundary is a terminus (TSS+), + strand transcript start. this orients the solver. this is critical! getting the orientation logic correct is a crucial
- boundary thus knows that the RIGHT sided message is composition-incompatible because *some* transcription ENDS at the boundary (some of the exon(1000,2000) transcription does not cross the boundary).
- boundary also knows that the LEFT RNA population DOES cross the boundary
- so the boundary composition equals the LEFT RNA abundance and the gDNA abundance 


We need to derive and design the logic and arithmetic for solving this time of exon|exon boundary. Let's call it "exon|exon boundary with terminus".

The key aspect of this is ORIENTING so we can figure out which transcript population is composition-compatible with the boundary. This requires:
- strand of the exons
- transcript start vs end

Do you agree?

Be sure that when you derive and design this, you teach me what you are doing, and ensure that I understand every piece of this. 


==================


Boundary at 2000:
- left region is exon (1000,2000)
- right region is exon+intron (2000,10000)
- the boundary is a splice junction (+ stranded)
- if message propagation is FORWARD (left -> right) then RNA SPLICES OUT (splice-out == leaving) at boundary 2000
- if message propagation is REVERSE (right -> left) then RNA SPLICED IN (splice-in == joining) at boundary 2000

So the logic depends on the direction of message propagation, the strand of the splice junction.






We need to augment our test chromosome to include more complicated cases with exon|exon boundaries. These boundaries start to emerge when transcripts have multiple isoforms.

the single-stranded exon|exon boundary is the first case to solve.

the behavior of the boundary depends on the following:
- splice junction? (yes/no)
- terminus? (yes/no)
- direction of message propagation (forward or backward)
- active strand (+ or -)

we need to be meticulous about the LOGIC at these boundaries.

here are examples:

- TA+ exons (2000, 3000), (10000, 13000)
- TB+ exons (2000, 3000), (11000, 13000)
- TC+ exons (12000, 13000)

At position 2000 we have an exon|exon boundary.




================================

# exon region -> exon|intron boundary message

Currently, we nullify these messages. Messages that are sent from exon region to exon|intron boundary (with a splice junction on the correct strand) must handle RNA fragments that SPLICE OUT.

The concept of "SPLICE OUT" implies that fragments do not cross the boundary as unspliced. Instead, they LEAVE via a splice junction. In many ways, this is the inverse operation from the "SPLICE IN" operation we have already solved.

From the standpoint of the exon REGION sending the message, there is nothing to do. The exon sends the message. Done.

We approach the exon region -> exon|intron boundary message from the perspective of the boundary that receives the message. One of the messages that the exon|intron boundary receives is from the EXON. That message contains purely unspliced fragments (RNA+, RNA-, gDNA).

The boundary measures spliced fragments and unspliced fragments. We'll proceed with this example with everything on the positive (+) strand for clarity.

The main question is, "how do we perform a composition transfer across the exon -> exon|intron boundary?"

We can compute the enrichment ratio as the ratio of total abundance:

`enrichment_ratio = exon_abundance / boundary_abundance`

Where `boundary_abundance = (spliced fragment abundance) + (unspliced fragment abundance). This is the same enrichment ratio that we computed for the SPLICE IN composition transfer.

Remember that `spliced fragment abundance` must be computed from all of the splice junctions that splice out. Same as what we did for splicing in.

First, we *rescale* the boundary by the enrichment ratio. This places the exon region and the exon|intron boundary on the same total abundance scale. Now, we can perform the transfer.

After the rescale, the spliced-out fragment abundance is subtracted from the exon's RNA abundance (on the same strand as the SJ). The subtraction is licensed because we have rescaled to the same total abundance.

After the subtraction, the leftover abundance is unspliced RNA + gDNA. That leftover abundance is still in on the exon's scale. We must rescale back to the boundary's frame of reference and the remaining abundance becomes the boundary's unspliced RNA + gDNA.

This is the arithmetic operation. This is an imputation. We also need to address the honest precision. 

The factors that affect precision are: 1) the exon count and effective length (abundance in counts/bp), 2) the boundary spliced fragment count and its length, 3) the boundary unspliced fragment count and length. Then there is the "cost" of imputation that we need to price as we have with other messages.
