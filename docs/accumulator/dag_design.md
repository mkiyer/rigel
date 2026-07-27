

# current behavior

currently rigel models the transcriptome as a linear graph of regions (genomic intervals) separated by boundaries (single positions).

The linear representation is simple and conducive to modeling fragments in contiguous genomic space. This representation allows straightforward modeling of genomic DNA (contiguous in genomic space) and unspliced RNA (contiguous in genomic space)

Multi-exonic transcripts have exons and introns. In the mature transcript, the introns are spliced out, leave the exons. The fragments crossing between exons are spliced fragments (pure RNA).

To accommodate spliced fragments, the boundary nodes have a compartment to store and track spliced mass.

The linear graph is bipartite: region <-> boundary <-> region <-> boundary. Traveling between two regions means visting a boundary.


## drawbacks of current model

- spliced fragments must be "split" across two boundaries. they contribute mass to two different boundaries
- this splitting conserves total fragment mass, but compromises fragment length. each splice junction boundary node has a very short fragment length
- fragments that span more than one boundary are fractured and partitioned across multiple nodes
- the discrete nature of counting is compromised.


# new behavior

- in the proposed design (this document), we wish to refactor our accumulator to store an acyclic graph representation. the linear (contiguous) graph is a special case of the acyclic graph.

- graph retains its bipartite nature. regions are genomic intervals. boundaries are transitions between regions. the architecture has not been determined yet.

- the previous linear model becomes ONE path through the graph. specifically, the graph retains the linear traversal in contiguous genomic space.

- splice junctions are fully represented as paths through the graph. they are no longer stored as separate information within a boundary node. they are separate paths through the graph.

- we wish to implement acyclic graph representation as a replacement for the current linear graph accumulator


# building the acyclic graph

the acyclic graph is constructed at index build time. it becomes part of the rigel index itself.

goal is to partition the transcriptome into 'regions' and 'boundaries'.

## region

- contiguous genomic interval e.g. (chr1, 1000, 2000)

## boundary

- single genomic position
- transcript start site (TSS)
- transcript stop site (TES)
- splice junction (exon-intron boundary)


## construction

construction reuses much of the same region partitioning code:
- loop over each chromosome
- sort transcripts by start, end positions
- enumerate all boundaries (TSS/TES, SJ positions)
- enumerate all splice junctions (introns)

### linear (contiguous genomic space)

- walk consecutive boundary pairs
- intervals between two boundaries are the regions
- edge between adjacent regions in genomic space is the contiguous genomic path (unspliced fragments)

### nonlinear (splice junctions)

- for each unique splice junction (intron), add an edge between two regions.
- **splice junctions are stranded (directed)** this can be represented as an edge direction, OR stored as an edge attribute. splice junctions have a genomic splice motif that defines the strand. the direction of the edge is the strand of the transcript.
- whether we model the graph as directed or undirected is an engineering design, not a modeling decision.

- Edges between regions in contiguous genomic space (linear) are undirected
- Splice junctions are stranded, therefore unidirectional


### nodes and edges?

with the new acyclic graph, it might make the most sense for 'regions' to be NODES in the graph, and 'boundaries' to be EDGES in the graph.

**note: we need to derive the most efficient graph representation**

### adding edges

The 'construction' above defines regions (nodes) and boundaries (edges) for each chromosome, but does not construct the edges.

- adjacent regions are always connected by boundaries. this captures the linear, contiguous representation that the current code uses. these are 'contiguous' boundaries. they are UNSPLICED boundaries.

- splice junctions are the new addition. splice junctions are *nonlinear* edges between boundaries


# new accumulator

the new acyclic graph representation is constructed by the index build. it should be stored/persisted as part of the rigel index. it should be loaded at runtime (rigel quant).

## software architecture comments

we need an efficient graph representation, likely in native code (C++). ideally, we would avoid a heavyweight dependency, unless there is an extremely straightforward, easy graph library that we can import. native speed is essential. the tool cannot slow down because of a transition to an acyclic graph accumulator.

## how to fragments translate to the graph?

we will need to rewrite our accumulator. there are some KEY PRINCIPLES:

### A FRAGMENT IS A PATH THROUGH THE GRAPH

previously, a fragment got broken into pieces if it crossed more than one regions, and its mass was stored inside the boundary nodes rather than the region nodes.

in the new model, we preserve fragments as atomic units that cannot be broken apart. the discrete unit of counting is the fragment. there are no fractional fragments.

this requires modeling a fragment as a PATH through the acyclic graph.

## construction of a path graph

- fragments are processed from the BAM file in the same way
- when fragments are DEPOSITED to the accumulator, our new accumulator takes over.

### algorithm for a single fragment

- map the fragment to a PATH in the graph
- the fragment has a strand, aligned blocks, and splice junctions
- this needs to be mapped to nodes (regions) of the acyclic graph, connected by a path through the graph.

#### reconciling 'implicit' splice junctions

some fragment alignments do not harbor the explicit splice junctions. the fragment can be mapped to a set of region nodes, but the splicing pattern may be missing from the fragment alignments. if there are multiple splicing combinations, the exact path of the fragment through the acyclic graph is not known.

the best way to reconcile fragments with implicit splicing patterns is to wait until we can approximate the rna fragment length distribution. then we can take all the compatible paths for that fragments and assign the fragment to the most likely path based on the fragment length distribution. for example, if some of the paths through the graph are >1000bp, they are unlikely to be possible given an observed fragment length distribution of 200bp.

this means that spliced implicit fragments with multiple possible paths would need to be held aside until RNA FL distribution was known, and then accumulated probabilistically. sort of like a "multimapping" fragment.

#### add simple multimapping support?

the current accumulator throws away multimapping reads. if we introduce a "two-pass" accumulator, we could consider recovering multi-mapping reads using a similar strategy. each multi-mapping fragment has multiple possible paths through the acyclic graph. we have no idea which is more likely during the first pass.

in a second pass, we could probabistically assign a multi-mapping fragment to one of the possible paths, using the counts of all the matching paths. same way we would handle a uniquely-aligned fragment that has multiple compatible paths.


## the path graph subsumes the splice graph

the splice graph built by the index does not store paths. it stores single edges between regions, where an edge can be contiguous in genomic space, or a splice junction to a nonlinear distant location.

paths sit on top of the splice graph. so how we add path information to the graph?


### worked example

worked example with transcripts on the positive + strand:

TA+ exons (1000, 2000), (10000, 12000), (20000, 22000)
TB+ exons (1500, 2000), (9000, 12000), (20000, 25000)
TC+ exons (1000, 2000), (5000, 5050), (6000, 6050), (10000, 12000), (20000, 25000)

The transcript 'TC' has two tiny exons (50bp each). A single fragment will not likely fit inside the region.

Regions:
- AA (0, 1000)
- A (1000, 1500)
- B (1500, 2000)
- C (2000, 5000)
- D (5000, 5050)
- E (5050, 6000)
- F (6000, 6050)
- G (6050, 9000)
- H (9000, 10000)
- I (10000, 12000)
- J (12000, 20000)
- K (20000, 22000)
- L (22000, 25000)
- ZZ (25000, end of chromosome)



Let's map fragments to paths:

1) Fragment: (1900,2000), (10000, 10100) spliced
- this fragment is not CONTAINED within a single region
- it maps to (B -> I). two regions, one edge
- so this fragment would accumulate as the path (B -> I), which is nonlinear (spliced) and positive stranded.

A fragment mapping to 2 regions could be stored as an EDGE between the two regions.

So the path (B -> I) would be accumulated as a single count on the edge between node B and node I.

2) Fragment: (1990, 2000), (5000,5050), (6000, 6050), (10000, 10090)
- crosses 4 regions, including the tiny exons
- mapping is (B -> D -> F -> I)
- landing zones are node B to node I.
- accumulate as an EDGE between node B and node I
- this is a DIFFERENT edge than the simple (B -> I). this is (B -> D -> F -> I).
- both are spliced. both are pure RNA. both start and end at the same "landing zone" regions.
- they imply different fragment lengths!
- the path (B -> D -> F -> I) contains region 'D' (50bp) and region F (50bp) plus the landing zones. So the effective length of this edge is different from the simple (B -> I) edge. This multi-splice edge has a longer fragment length 

3) Fragment (1800,2050)
- contiguous in genomic space, NOT spliced
- maps to (B -> C)
- accumulate as edge between B and C.


### what do the worked examples show?

- goal of calibration is to estimate gdna vs rna
- NODES (Regions): accumulate CONTAINED fragments
- EDGES (Boundaries): fragments that cross one or more region boundaries

### how does accumulator work?

- convert fragment -> path (set of nodes)
- if fragment maps to single node, add it to the node
- if fragment maps to multiple nodes, then it adds as an edge between the first and last node.


### how to we reconcile paths vs single nodes?

**NOTE: we need derivation work here**

if we have contiguous genomic space with small nodes, for example, due to multiple transcription start sites or stop sites, then we can have unspliced fragments that cross multiple linear nodes..

*unspliced fragment maps to multiple nodes in contiguous space*

ex. contiguous nodes A = 10bp, B = 50bp, C = 40bp, D = 1000bp

We could have a fragment mapping to (A -> B -> C -> D). Individual nodes are too small to 'contain' fragments. We could also have fragments (B -> C -> D), and (C -> D). These are subpaths of (A -> B -> C -> D).

How do we handle the unspliced fragments with different subpaths that share the same contiguous (linear) path? We cannot collapse (A -> B ->C -> D) into a single node because they may have spliced junctions, etc. In our current code, we store these as "boundary" nodes even though they are not spliced. So the individual regions have ZERO contained fragments. Does this still work?

If we accumulate PATHS, how do we project them onto the same plane for calibration? 


### how do we maintain count conservation

if a fragment overlaps (A -> B -> C -> D) in contiguous genomic space, we do not want to credit each of the regions with +1 count. in the current code, the fragment becomes a boundary-crossing fragment and credits (A->B), (B->C), and (C->D) each with 0.33 count mass. If we accumulate "paths", we have a single edge (A -> B -> C -> D).. this is an edge from node A -> node D, unspliced, +1. Is this the correct way to do the accumulation?


### can we collapse the paths?

if all the matters is the starting node and the ending node, then for a fragment mapping to nodes (A -> B -> C -> D), we would merely need to store:

- edge src A dst D unspliced +1 count

this avoids the needs to store the internal path information. the question is, how much does the internal path matter? what really matters is: "is this spliced or unspliced?".

i think this is probably correct!

but there is a catch -- the fragment lengths matter (for unspliced fragments).

# TSS and TES

- a major advantage of the new model is that we will maintain TSS and TES information.

- transcript start site (TSS) must be added as a SOURCE node to the splice graph
- transcript end sites (TES) becomes SINK node

- TSS and TES are the start/end of our forward-backward belief propagation at calibration time

- if calibration processes one locus at a time, a single locus can share a single TSS node with edges from the TSS node to all of the start site nodes (directional)

- a single locus can have a single TES node with edges to all the transcript end sites



# store fragment lengths?

we were reluctant to store a fragment length distribution at every unspliced node or edge. we don't need to store spliced fragment lengths (we know they are RNA already, and they are used to train the RNA FL distribution). similarly we don't need to store intergenic fragment lengths (they are pure gDNA and train the gDNA FL distribution).

But we could decide to store unspliced fragment lengths (as a histogram, even a heavily quantized histogram for memory). this gives us a pre-node deconvolution signal


# comments and questions as of first draft


- currently, any fragment that overlaps >1 region is accumulated as an edge between two regions. a fragment may overlap many regions, but is still accumulated as a single exon between its starting and ending regions. is this acceptable? are we discarding critical information?

- this accumulation question applies to both spliced and unspliced fragments, which can splice AND/OR cross regions contiguously. a single fragment could cross many regions, some via splicing, some via contiguous crossing. a single splice junction is all that is required to render the entire fragment 'spliced'.

- this model would be simple. it would give us what we already have, except more elegant, and give us the splicing pattern information as well.

- at a single node, we could collapse the edges into a 'spliced counts' and 'unspliced counts' entering/exiting the node. this is exactly what we have now, without the fractional mass allocation and mass-splitting behavior.

- what happens to 'tiny' nodes? small regions cannot harbor contained fragments. they will be sparse or zero, but they will remain connected to other nodes via the acyclic graph

- empty nodes will 'relay' messages during belief propagation


















