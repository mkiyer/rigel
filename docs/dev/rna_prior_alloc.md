# Overview

Prior to this session we had a single grouped RNA prior. The RNA prior is UNSPLICED fragments (not total fragments). The RNA prior is then applied to all transcripts. Was our initial strategy to apply the prior uniformly to all transcripts? We had the warm start based on coverage weights, and then the prior multiplies by the warm start. Is that correct?

Then, we added a hard ZERO to the synthetic nascent RNA. They can still get a warm start, but ZERO prior mass.

Now, we need to design a more intelligent method. Our grouped UNSPLICED RNA prior is obtained from calibration, after deconvoluting genomic DNA.

There are several issues that this session must address:

1) calibration has region/boundary granularity (each region/boundary has a gDNA vs RNA compositional breakdown). We have data that could estimate unspliced counts per transcript, but don't use it.

2) transcript GEOMETRY governs how many unspliced fragments it can support. If a transcript has many tiny exons, all of its counts could be spliced. If a transcript has a few large exons, all of its counts will be unspliced. We do not regard transcript geometry

-----

The goal: after calibration finishes, we run "assemble priors" per locus. Assemble priors produces a gDNA prior in pseudocounts, and an RNA prior in pseudocounts. The gDNA prior will not change.

In this new version, we will partition the RNA prior pseudocounts across the different transcripts. Each transcript will be given a portion of the RNA prior pseudocounts relative to its calibrated regions/boundaries.

The challenge is HOW to assign a per-transcript prior

----

# Transcript to REGION/BOUNDARY map

A transcript consists of exons (contiguous genomic intervals.We already have a method to map a transcript onto the splice graph (regions are genomic spans, boundary boundaries are points between regions, splice junctions are jumps between two boundaries).

Each transcript must be mapped to:
- regions (contiguous genomic spans completely overlapped by the transcript)
- splice junctions that are used by the transcript
- boundaries that the transcript crosses contiguous (on both sides)

This already exists in the `_transcript_node_incidence` function.


---

# Transcript upper bound density can be computed

If we transform a transcript to its region/boundary map, the transcript becomes a SEQUENCE of (region/boundary/sj) objects.

If we iterate along the sequence of region/boundary/sj objects:

1) REGIONS: have a contained mass equal to integer count of fragments contained
2) BOUNDARIES: have a fractional mass, if a fragment overlaps multiple boundaries its mass its 1.0 count shared among its boundaries
3) SPLICE JUNCTIONS: have a fractional mass, if a fragment crosses multiple splice junctions/boundaries its 1.0 count is portioned among the these sj + boundaries.

Remember, the mass along this sequence of region/boundary/sj is SHARED among all transcripts that are compatible. The deconvolution has been done for gDNA vs RNA, but not per-transcript yet.


# Theorem for transcript prior allocation

Our theorem is as follows:
- A transcript's upper bound density is the MINIMUM of densities of its region/boundary/sj objects.

However, a naive implementation of this theorem will fail catastrophically. This is because fragments can land sparsely along a transcript. Any region/boundary/sj with a low randomly sampled density will crush the entire transcript, just due to sampling variation. Sampling variation is worse when regions/boundaries/sj are small. Boundaries and SJ have a short length. So the transcript's entire density hangs in the balance of the weakest, lowest count boundary/sj along its entire path.


# Naive solution for transcript prior allocation

Our first prior allocation method aims to balance simplicity with accuracy.

We have several candidates:
- harmonic mean along transcript region/boundary/sj
- geometric mean
- arithmetic mean

Given the theorem above (under uniform coverage with no coverage undulation/variation, the transcript density is the minimum density along its path), we would prefer the HARMONIC MEAN.

However, the harmonic mean is undefined if there is a ZERO along the transcript's path. By our theorem, a ZERO along a transcripts path means it cannot be expressed, but the difference between a true zero and a false zero could simply be random sampling error. This is where we can raise precision as a possible factor. If we have a zero over 10000bp, that is a more confident zero than a zero over 200bp.

Open question:
- should the mass along the path be precision-weighted? So imprecise estimates are down-weighted relative to more precise estimates?

Options for the harmonic mean: 
- add a small constant to 'rescue' zeros
- omit regions/boundary/sj that are zero

Remember that this method approximates the MASS for each transcript. The MASS is different from density. 

Density is computed by an effective length computation. We use an 'effective length shrinkage' method that accounts for enrichment/depletion by hybrid capture. So regions/boundary/sj that are ZERO may be zero because they are off-target (no probes to capture them). 

The analogous approach for the harmonic mean would be to omit regions/boundary/sj that are ZERO. Our effective length method also ignores regions/boundary/sj that are zero.. they are not usable because we don't know whether it is a true zero or an off-target/depleted zero.

In summary:
- We compute a per-transcript MASS which is the harmonic (geometric?) mean over the SEQUENCE of regions/boundaries/sj in its path. This becomes a WEIGHT for each transcript.

**KEY REALIZATION**

- Our RNA prior is UNSPLICED fragments only, because spliced fragments are already handled by the EM (they don't compete with gDNA).

- Therefore, the transcript mass harmonic mean method must exclude the spliced mass. The weight must only include the unspliced mass.

- After each transcript has a weight, then we use the per-transcript weights to allocate the total RNA prior pseudocounts in a weighted fashion.


---

# splice junction index and accumulator

the rigel index already builds a splice junction index, where each splice junction is defined by:
- ref (reference chromosome)
- strand (splice junctions are stranded)
- start (genome start of the intron)
- end (genome end of the intron)

there is redundancy here, because biology does not support (+) and (-) splice junctions at the same exact point, in otherwords, the following is IMPOSSIBLE:
- chr1, +, 1000, 2000
- chr1, -, 1000, 2000

nevertheless, we store this for completeness and error checking, as we need all four fields to uniquely define a splice junction


## BAM scan / accumulator

Rigel scans a BAM file and builds Fragment objects. Fragments with splice junctions are checked against artifacts (a blacklist database). Fragments can also have "implicit" splice junctions in the unsequenced gap between to reads of the fragment. These can all be resolved during the BAM scan phase.

Fragments are then "deposited" to the accumulator. At the time of deposit, the spliced fragments are fixed.

The rigel accumulator should maintain a splice junction accumulator that tracks spliced fragments:

- ref, strand, start, end (the sj key)
- counts[2] (2-integer vector with counts on each aligned strand)
- mass[2] (2-float vector, a spliced fragment may cross multiple sj but its mass must be conserved and shared, so that its 1.0 count is divided across all of the splice junctions the fragment crosses. So the mass of a spliced fragment equals 1.0 if it only crosses one sj, otherwise, the mass is <1). We allocate the mass proportionally based on the overlap bases related to either sj. Spliced sj may cross OTHER boundaries (contiguously, not spliced), which can also decrease their overall mass. This guarantees conservation of fragment mass.

---

# Boundaries versus Splice Junctions

The true theoretical rigel index is a graph. Regions are contiguous genomic intervals. Boundaries connect regions.

The current implementation uses shortcuts and avoids building a formal graph structure. Instead it uses genomic contiguity to index regions and boundaries.

So for example:
- region A (chr1, 1000,2000)
- region B (chr1, 2000,4000)
- region C (chr1, 4000,6000)
- region D (chr1, 6000,10000)

Regions are contiguous across the genome. There is ALWAYS a contiguous boundary between two regions: for example, there is boundary at position 2000 (A <-> B), 4000 (B <-> C), and 6000 (C <-> D). Rigel shortcuts building these boundary structures because the genome is contiguous the next index is the contiguous region.

So, currently rigel refers to 'boundaries' as boundaries between two contiguous regions.

However, splice junctions are also 'boundaries', but they are between discontiguous regions. 

For example, we could have a splice junction:
- sj AC+ (chr1, +, 2000, 4000)

This becomes an unidirectional boundary connecting A -> C. 

---

In order for us to compute our per-transcript priors, we need to access the individual splice junction data for each transcript.

A transcript will map to an explicit sequence of regions and boundaries. We need to be able to lookup these counts EXACTLY.

We need a detailed architectural audit to ensure we have the efficient data model and accumulator structures to support this project.

What data structures do we need to modify and/or build?

---

# the splice junction counts and mass as a primary tool output

some investigators want to access the raw spliced counts and spliced mass associated with the data. this is sometimes considered a version of the 'true' per-transcript counts, because spliced fragments typically can be trusted to be pure RNA. this should be a primary output of the rigel tool.

so building this correctly in the accumulator phase is worthwhile and valid.

---



---



## IGNORE THE BELOW -- this is a much more complex idea for a smoothing approach following by transcript enumeration. This is incomplete

** IGNORE THE BELOW **


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

In order to do the smoothing, every region/boundary/sj needs to track the ends of fragments, so that the algorithm can extend fragment ends upstream (5', toward the TSS) and downstream (3' toward the TES).

**The accumulator does not currently track fragment ends** It would need to store 'count_left' and 'count_right'. 









