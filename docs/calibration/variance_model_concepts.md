
# Variance / precision model concepts

RNA-seq is a mixture of genomic DNA (gDNA) and RNA. Rigel deconvolutes ambiguous fragments into gDNA and RNA pools.

The total unspliced counts within a node are observed and fixed. The solver must assign the total unspliced count budget to DNA and RNA.

gDNA and RNA have different fragment length distributions. Therefore, ABUNDANCE (also called density) of DNA and RNA in a node depends on their relative composition. For example, given a count budget of 1000 unspliced counts, DNA FL mean 250, RNA FL mean 100, the ABUNDANCE (density) of the node can vary between 1000/250 = 4 (pure DNA) and 1000/100 = 10 (pure RNA).

The fragment length distributions are modeled and known (earlier step of calibration).

## unified solver pathway

We are working on the new solver, called the "unified solver". We are trying to retire and delete the old code paths. Focus all work on the new solver.


## Initialization

At the start of calibration, we setup any nodes that are known (measured counts), and then derive nodes that can be solved (intron deconvolution, strand deconvolution).

### MEASURED counts get poisson precision

- intergenic nodes are pure gDNA, known
- intergenic-exon nodes are pure gDNA, known
- spliced fragments in boundaries are pure RNA, known

There are two sources of directly measured counts:
- intergenic fragments are pure gDNA
- spliced fragments are pure RNA

The SOURCES of measured biology get count precision:
- Poisson now (simple)
- Negative binomial later (deferred), use simulator generates poisson counts. we have not derived a method to fit overdispersion yet

### "intron factory" - deriving gDNA from introns

- Fragments in intron nodes = gDNA + unspliced RNA
- These intronic fragment can be deconvoluted based on the intergenic (pure gDNA) distribution
- We have a negative binomial model for this
- intron nodes estimate their gDNA content using this approach. this comes with a precision based on our NB model.

- the gDNA within introns can be deconvolved by comparing to intergenic distribution, without strand information. 

- However, the residual RNA within introns requires strand assignment. if the intron is single-stranded, strand can be inferred. if the intron is ambiguous (introns of two transcripts on opposite strands overlap), then the RNA strand cannot be inferred. the ambiguous/both-strand case can STILL peel out the gDNA, but cannot assign strand to the RNA.


### strand-specific deconvolution

- nodes with 1-DOF can directly solve with strand specific data alone using our strand balance deconvolution model
- solved nodes get precision set by the strand model


### unsolved nodes - set zero precision defaults

- many nodes cannot be directly solved. whether or not solve is possible depends on whether or not we have strand-specific data or unstranded data.

- "measured" nodes are 0-DOF and always solved



#### strand-specific data

- single-stranded nodes can be directly solved with strand-specific data
- both-stranded (ambiguous) nodes cannot be fully solved with strand-specific data. However, they CAN be partially solved!
- the strand balance provides what we call the "tilt" of a node. This gives us 1-DOF. Ambiguous nodes have 2-DOF which we model as "tilt" (one DOF) and "gdna content" (2nd DOF). Knowing the strand tilt CONSTRAINS the solution to two possibilities. It bounds the possible range of gDNA.




### Summary of initialization

Initialization sets up ndoes wtih 


