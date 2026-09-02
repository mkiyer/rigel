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









