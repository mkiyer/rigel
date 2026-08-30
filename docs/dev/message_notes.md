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


