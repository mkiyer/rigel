# Overview

Prior to this session we had a single grouped RNA prior. The RNA prior is UNSPLICED fragments (not total fragments). The RNA prior is then applied to all transcripts. Was our initial strategy to apply the prior uniformly to all transcripts? We had the warm start based on coverage weights, and then the prior multiplies by the warm start. Is that correct?

Then, we added a hard ZERO to the synthetic nascent RNA. They can still get a warm start, but ZERO prior mass.

Now, we need to design a more intelligent method. Our grouped UNSPLICED RNA prior is obtained from calibration, after deconvoluting genomic DNA.

There are several issues that this session must address:

1) calibration has node/edge granularity (each node/edge has a gDNA vs RNA compositional breakdown). We have data that could estimate unspliced counts per transcript, but don't use it.

2) transcript GEOMETRY governs how many unspliced fragments it can support. If a transcript has many tiny exons, all of its counts could be spliced. If a transcript has a few large exons, all of its counts will be unspliced. We do not regard transcript geometry

-----

The goal: after calibration finishes, we run "assemble priors" per locus. Assemble priors produces a gDNA prior in pseudocounts, and an RNA prior in pseudocounts. The gDNA prior will not change.

In this new version, we will partition the RNA prior pseudocounts across the different transcripts. Each transcript will be given a portion of the RNA prior pseudocounts relative to its calibrated nodes/edges.

The challenge is HOW to assign a per-transcript prior

----

# Transcript to NODE/EDGE map

A transcript consists of exons (contiguous genomic intervals.We already have a method to map a transcript onto the splice graph (nodes are genomic spans, boundary edges are points between nodes, splice junctions are jumps between two edges).

Each transcript must be mapped to:
- nodes (contiguous genomic spans completely overlapped by the transcript)
- splice junctions that are used by the transcript
- edges that the transcript crosses contiguous (on both sides)

This already exists in the `_transcript_node_incidence` function.

---

# Transcript upper bound density can be computed

If we transform a transcript to its node/edge map, the transcript becomes a SEQUENCE of (node/edge/sj) objects.

If we iterate along the sequence of node/edge/sj objects:

1) NODES: have a contained mass equal to integer count of fragments contained
2) EDGES: have a fractional mass, if a fragment overlaps multiple edges its mass its 1.0 count shared among its edges
3) SPLICE JUNCTIONS: have a fractional mass, if a fragment crosses multiple splice junctions/edges its 1.0 count is portioned among the these sj + edges.

Remember, the mass along this sequence of node/edge/sj is SHARED among all transcripts that are compatible. The deconvolution has been done for gDNA vs RNA, but not per-transcript yet.


# Theorem for transcript prior allocation

Our theorem is as follows:
- A transcript's upper bound density is the MINIMUM of densities of its node/edge/sj objects.

However, a naive implementation of this theorem will fail catastrophically. This is because fragments can land sparsely along a transcript. Any node/edge/sj with a low randomly sampled density will crush the entire transcript, just due to sampling variation. Sampling variation is worse when nodes/edges/sj are small. Edges and SJ have a short length. So the transcript's entire density hangs in the balance of the weakest, lowest count edge/sj along its entire path.


## A two stage algorithm

### Stage 1: smoothing

- Assume every fragment is sampled from a single transcript (or gdna)
- The fragment does not span the entire transcript/gdna because of the sequencing platform (short reads)
- In order to apply our 'minimum along transcript path' method, we must extend each fragment in both directions until it reaches the start and end of a transcript.
- A single fragment may be compatible with MANY transcripts, and extending the fragment must be done probabilistically.

Here's an example:

TA+ exons (1000, 2000), (9000, 10000)
TB+ exons (500, 11000)
TC+ exons (500, 2000), (5000, 6000), (9000, 11000)
TD+ exons (4000, 6000), (9000, 11000)

Take a 200bp fragment at (5500, 5700). It is compatible with TB+ and TC+. The fragment must be extended to the left (lower genomic coordinates) and to the right (higher genomic coordinates) until it reaches TSS (transcript start) and TES (transcript end).

In order to do the smoothing, every node/edge/sj needs to track the ends of fragments, so that the algorithm can extend fragment ends upstream (5', toward the TSS) and downstream (3' toward the TES).

**The accumulator does not currently track fragment ends** It would need to store 'count_left' and 'count_right'. 









