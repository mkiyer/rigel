
# New enrichment-weighted effective length method

## Setup:
- Locus genomic interval (ref, locus_start, locus_end)
- Transcripts (both mature and nascent RNA)
- Mature RNA defined over its EXONS (introns excluded)
- Nascent RNA defined over its full span

## Genomic effective length: every component has a genome-defined effective length (sometimes referred to as the geometric effective length).
- Mature RNA: transcript length (exons only) minus integral of RNA FL distribution (clipped at the transcript length). There should be code for this.
- Nascent RNA: transcript genomic span minus integral of RNA FL distribution (clipped at total transcript genomic span). Same code for mature RNA
- gDNA: locus genomic span (we do not subtract the gDNA FL distribution here because gDNA extends beyond the locus)

## Seed regions

- Intergenic regions are considered pure gDNA. There are relatively rare unannotated RNAs that can bias intergenic gDNA in local areas, but intergenic regions globally should be approximately 'pure' gDNA
- Nascent RNA is present in introns, but tends to be very sparse. gDNA is ~uniform. Therefore, intronic regions should roughly approximately gDNA
- Intergenic <-> Exon boundaries can contain gDNA that is captured by hybrid capture probes, and are therefore an important entity and a type of 'seed'
- Intronic <-> Exon boundaries are often captured by hybrid capture probe panels and similarly estimate gDNA. Intron <-> Exon boundaries include (sparse) nascent RNA. Globally, we assume that Intron <-> Exon boundary fragments approximately gDNA


## Global gDNA density

Global density should start as a global scalar, derived from any observable seed regions. This is the sum of the gDNA mass divided by the cumulative length (bp) of the regions.

How do handle hybrid capture and the global gDNA estimate? This must be data driven. We build nonparametric var ~ mean models that work well with hybrid capture. 

At initialization, we only have our seed regions/boundaries. A subset of the intergenic <-> exon and intron <-> exon boundaries could estimate captured gDNA. The problem is that probe panels typically target a fraction of the transcripts. Sometimes a very small fraction of transcripts are targeted. Therefore we cannot assume that all EXONS are captured, or our model will be built incorrectly. Our local estimates uses boundary crossing fragments. Our initial global estimate should be 'everything is uniform'. We may need to iterate so that the global gDNA estimate reflects ALL the regions/boundaries rather than just the biased seed nodes.

We can consider *ITERATING* -- after calibration finishes, we have a gDNA estimation at all regions/boundaries, not just seed regions. We could re-initialize the global gDNA prior and rerun.


## Nodes: Regions and Boundaries

There are two node types:

- Regions: genomic intervals (ref, start, end) that store 'contained' fragment mass

Region geometric effective length is therefore constrained by the region size and the fragment length distribution.

```
region_geom_eff_len = region_len - frag_len
```

Where 'frag_len' is the integral of the FL distribution (gDNA or RNA) clipped at the region_size (fragments must be shorter than the region size to be contained)


- Boundaries: single genomic positions (ref, pos) that store crossing fragment mass (flux).

Boundary geometric effective length is the fragment length distribution (boundaries can be thought of as regions of length 1)

```
boundary_geom_eff_len = frag_len - 1
```

Where fragment length is the integral of the FL distribution (gDNA or RNA)

Boundary density must be computed carefully. The regions on each side of the boundary constrain the boundary mass (fragments can cross multiple regions and multiple boundaries and we enforce conservation of fragment mass).

Functions to compute these values should already exist.

**KEY INSIGHT**

**Boundaries do more than measure flux, they own mass**

Boundary crossing fragments are NOT stored in 'regions'. Fragments that crossing boundaries are 'owned' by the boundaries. *This makes boundary nodes first-class fragment mass carriers*. . Their mass is divided into two 'sides'. Each 'side' contains a portion of the total mass based on base overlap.


## Enrichment ratios

The global gDNA mass devided by the global genome size (bp) is the global gDNA density.

An `enrichment ratio` is defined as the ratio of local gDNA density (node/transcript/locus) to global gDNA density.
- If Enrichment ratio > 1.0: the node/transcript/locus is enriched
- If Enrichment ratio < 1.0: depleted


### Computing enrichment ratios

- Single node (region/boundary):  Computation of enrichment ratios for single nodes (regions or boundaries) is straightforward. We divide the node gDNA density by the global gDNA density.

- Component (locus, transcript): This requires a component -> node map. If the component is a transcript, we must map the transcript onto a collection of regions and boundaries.

Given transcript TX:
- Find leftmost region start (note there was a major bug here in the past, ensure the leftmost start is correct)
- Find the rightmost region end (same concern, ensure correct, regression test)

- Collect regions and boundaries. Note that boundaries have TWO SIDES. 

Transcripts must use the appropriate boundary sides.
    - For a transcript, the first/last regions should not use the END boundary (boundary crossing fragments that overhang the transcript are not compatible with the transcript span). This includes the LEFTMOST boundary of the LEFT region and the RIGHTMOST boundary of the RIGHT region
    - For a transcript, exon <-> intron boundaries must be handled carefully. We include only EXON regions (use region signature) and exclude INTRON regions. For exon <-> intron boundaries, we only use ONE SIDE of the boundary (the side touching the exon). The other side of the boundary touches the intron.

We need to CAREFULLY build the transcript -> node map
    - ensure we map a transcript to the appropriate nodes (regions and boundaries) and the appropriate boundary exons (respecting introns)
    - regression tests will be needed
    - careful here! no bugs

Start/end overlap fix:
- Regions are not guaranteed to EXACTLY overlap a transcript span. Multiple transcripts on the same strand with slightly different start sites can be grouped into a single region.
- We could ADJUST a transcript enrichment ratio calculation to account for differences in region overlap.

Once a transcript is mapped onto nodes (regions/boundaries), an enrichment ratio can be computed:
- sum of gDNA mass across nodes (numerator)
- sum of gDNA lengths across nodes (denominator)




