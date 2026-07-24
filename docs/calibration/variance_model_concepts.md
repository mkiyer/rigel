
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

- single-stranded nodes can be directly solved with strand-specific data
- both-stranded (ambiguous) nodes cannot be fully solved with strand-specific data. However, they CAN be partially solved!
- the strand balance provides what we call the "tilt" of a node. This gives us 1-DOF. Ambiguous nodes have 2-DOF which we model as "tilt" (one DOF) and "gdna content" (2nd DOF). Knowing the strand tilt CONSTRAINS the solution to two possibilities. It bounds the possible range of gDNA.


### unsolved nodes - set zero precision defaults

- nodes may or may not be "solvable" at initialization time.

- solvability depends on 1) node structure (which RNA strands are active). Intergenic = 0-DOF. Single-stranded = 1-DOF. Both-stranded == ambiguous = 2-DOF. 
- Intronic nodes (absence of exons) is a special case that allows us to determine gDNA by comparing directly to the intergenic distribution. Single-strand introns allow estimation of gDNA AND RNA (inferred). Both-stranded introns (2-DOF) can still estimate gDNA but cannot assign strand to RNA. Strand-specific data allows a FULL SOLVE of 2-DOF introns.

- Nodes MAY NOT be solvable at initialization.
- Nodes that are not solved must be set to default values.
- The default values are 100% gDNA, ZERO precision.


### Summary of initialization

Initialization sets up calibration to perform a "prior-free" solve. We call this "pass0". The prior-free solve estimates initial values nodes through forward-backward belief propagation (message passing).

Our goal is to formalize, solidify, and harden the initialization phase in our unified solver.

- WE need unit tests for each piece described above.
- We need to know that initialization works as described
- The code needs to be completely cleaned up and made clear, concise, readable, and maintainable.
- Once initialization is hardened and ready, the next step is to derive the variance model for the pass0 solve.


## Variance model

Our variance model needs to be rebuilt from the ground up:

1) define sources of variance
2) model sources of variance

Goals of the variance model are:
1) avoid assumptions about the distribution of the data. we are modeling simulated data, but real data has complex sources of variation
2) honest precision based on discrete counts measured over a discrete length (base pairs)
3) precision != accuracy. High counts are more precise, but this does not make the DECONVOLUTION result more accurate. This is the "count-zero" precision paradigm. Only measured-pure counts get count precision


### MEASURED counts get count precision

- intergenic fragments -> pure gDNA -> count precision (poisson now, NB later)
- spliced fragment -> pure RNA -> count precision (poisson now, NB later)


### Strand deconvolution precision

- our strand model is modeled using beta binomial. overdispersiom is fit to the data.
- nodes solved with the strand solver get count precision

** NOTE ** we found that strand balance measurement accuracy dramatically impacts solver performance. With unstranded data, our rna_sense_fract ~ 0.50, but in reality it is never exactly 0.50. tiny deviation from 0.50 causes cascading phantom precision problems. We have SOLVED this. This is a key aspect of our strand solver that needs to be incorporated.

### Density deconvolution precision

In the absence of strand-specific data, we can deconvolute regions using abundance (counts/length, or density).

However, density deconvolution requires a gDNA prior. IF we know the distribution of gDNA, we can model unspliced counts as:
- total density = rna density + gdna density 
- if gdna density can be estimated, then total density can be estimated.

we must respect fragment length here!

total density depends on composition of rna + density:
- rna density = rna counts / rna FL
- gdna density = gdna counts / dna FL
- total density = rna_density + gdna_density

We need to build a generic count deconvolution system. fortunately, we already have one called the "intron factory", but its disguised as a special case of the generic density deconvolution problem.

The gDNA prior (or hyperprior) acts as the floor for the solve. The residual RNA must be allocated to both strands. single-stranded nodes have 1-DOF and so the RNA can be inferred. both-stranded (AMBIG) nodes have 1-DOF and requires additional information to fully solve (such as strand-specific data, or message propgation).


### Intron deconvolution precision

- Intron deconvolution is a special case of density deconvolution where we KNOW the gDNA hyperprior at initialization time 

- This requires the assumption that **INTERGENIC NODES AND INTRONS ARE NOT CAPTURED BY HYBRID CAPTURE PROBES**. If this assumption is violated, then we CANNOT solve introns at initialization time.

- The assumption is violated if we ever have a hybrid capture probe panel that places probes in introns or intergenic regions (this is rare for most real rna-seq data).

- With this assumption, the gDNA hyperprior simplifies to the intergenic node distribution. The intergenic nodes measured the gDNA background level (depleted under hybrid capture, uniform without hybrid capture)

- Therefore, we can solve introns with a default gDNA hyperprior: the intergenic node distribution.

- Introns are deconvoluted into gDNA + unstranded RNA using a negative binomial model. This assigns honest precision to the deconvoluted components


### Combined intron deconvolution + strand deconvolution -> unified solve

There is SYNERGY between the strand deconvolution and the intron deconvolution approaches. We should NOT think of them as independent solving strategies.

Together, intron deconv + strand deconv can directly solving a 2-DOF (intron) node.

- We need to properly handle the unknown RNA strand. When the intron is single-stranded, the node is fully solved (1-DOF).

- When the intron is both-stranded (ambig), we cannot assign strand to the RNA EXCEPT if have strand-specific data. Remember, there are 2-DOF for ambiguous nodes -- Strand-specific data provides the "tilt", and the intron deconvolution provides the gDNA level, so these nodes are solvable.

- We need an ELEGANT integration of strand solve + intron deconvolution that fully solves ambiguous introns and assigns them honest precision (precision of strand deconvolution, precision of intronic gDNA deconvolution)

- this needs to validated and tested extensively. with regression tests.

- we need unified intron solve with merged precision of the strand deconv and intron deconv approaches


## "Composition" precision

### What is a "solve"? What is the solution?

Each node must partition its UNSPLICED counts across the active components (RNA+, RNA-, gDNA). 

A "solve" consists of:
1) a partitioning of unspliced counts to (RNA+, RNA-, gDNA).
2) precisions associated with the partitioning

Total unspliced counts are fixed:
- total unspliced counts = rna+ counts + rna- counts + gDNA counts.
- therefore, the solution can be thought of as a pie chart that we are allocating to three pie pieces.

We only need 2-DOF to model the solution, because the three components share the same total unspliced count pool.

### Solve state:

1) the strand tilt (strand balance)
2) the gDNA level

Two fractions [0,1] each with a precision [0, infinity].

**OPEN QUESTION** Do you agree with this derivation? Is this what our solver currently models?

### Goal of initialization

Each node must have a "current belief" that includes its state (two fractions, two precisions).

Once this framework is hardened, we can begin to model the variance of message propagation.

========================

You next task is to review survey the existing variance source. 

WE need to clean up the code, consolidate the different variance pieces, and ensure that the code is clean, concise, readable, maintainable.

1) Ensure clean models for the sources of variance (measured counts, strand deconvolution, density deconvolution)
2) Ensure we correctly model the SYNERGY between density deconvolution and strand deconvolution, a form of "self-solve"

WE need to consolidate existing code, clean it up, and synthesize everything into a clean, clear format.

We must be able to "solve" all of the above cases. None of these involve message propagation. That comes later. For now, we have:
- measured density with count precision
- strand deconvolution (with strand specific data) with BB precision -> ENSURE we handle the unstranded cases properly (this was previously handled well)
- density deconvolution (with a gDNA prior) with negative binomial precision (intron deconvolution is the current case that we use this)

- integrated / synergistic strand+density deconvolution given strand-specific data and gdna hyperprior.
-- the integrated solver is our goal!
-- handles the full spectrum of strand (from unstranded to perfectly strand specific)
-- handles case where gDNA prior is available (introns) or unavailable (exons). AFTER pass-0 solve, we will construct the gDNA hyperprior from the solved nodes, enabling integrated deconvolution of previously unsolvable nodes.

**WE DO NOT HANDLE MESSAGE PROPAGATION YET**
**THIS IS DEFERRED, BECAUSE MESSAGE PROPAGATION REQUIRES MERGING/COMBINING PRECISIONS FROM VARIOUS COMPONENTS**

=====


