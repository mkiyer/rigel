# TODO


## Gap SJ correction

I would like to review our "Gap splice junction" (Gap SJ) algorithm and concept. There might be an opportunity for improvement. Gap splice jucntion detection and subtraction is intended to solve cases where a fragment is not completely sequenced by paired-end reads leaving an unsequenced "gap" in the center of the fragment. If the gap contains splice junctions (introns) the splice junctions will not be detected and the fragment may appear to be unspliced (with a very long fragment length). Our current approach looks for introns within the gap that can be subtracted from the fragment length, effectively correcting the fragment. For the purposes of scoring and competition in the EM, these fragment I believe still remain "unspliced". So we rely on the fragment length calculation to allow the mature RNA candidates to outcompete nascent RNA or genomic DNA in the EM. IF we correct for Gap introns/SJ correctly, the mature RNA candidates should score MUCH MUCH better than their nascent RNA or gDNA equivalents. Now, I have another idea for gap SJ correction. At the point where we are resolving gap splice junctions, we have both a fragment and its transcript match. We KNOW the transcript structure. So if we translate the fragment genomic coordinates to transcript coordinates, we could then compute the fragment length trivially in transcript coordinate space (rather than genomic coordinate space). We would find the transcript coordinates of the fragment and compute the fragment length in transcript coordinate space. The rigel tool should already have efficient code to convert from genomic to transcript coordinates -- this could be reused or we could develop an extremely computationally efficient method for this. If we compute fragment length this way, then the fragment length computation will be robust to tiny overhang issues (where we few bases overhang into intronic space). --- I would like you to evaluate this carefully and propose a plan. Do you think this is a better idea than our current implementation? Will it be more robust? If so, create an implementation plan to replace the existing gap splice junction detection and subtraction code. I'll review the plan before we implement it.





## Unannotated splice junctions

A fraction of splice junctions are "unannotated" in that they don't have exact matches to the reference. These can be either artifacts (aligner errors) or 'novel' isoforms/transcripts. There is a lot of information that we are not taking advantage of yet:

- partially annotated: 3' site, 5' site, or both
- (easy) is the 5' splice site annotated
- (easy) is the 3' splice site annotated
- (medium) how much 'anchor' on either side of the intron?
- (hard) does this junction match any of the 'blacklist' splice junctions -- requires building a splice junction blacklist.

### False positive splice junctions

We can "blacklist" splice junction artifacts by doing genome-wide mapping. Align to the genome and detect cases where gDNA sequence -> spliced alignment (artifact). Compile a "blacklist" from these cases.



## Evidence-guided Nascent RNA gating

- We could require evidence and gate nascent RNAs
- Evidence would be fragment within an 'unambiguously intronic' genomic interval
- In the new model, nascent RNAs are just (synthetic) transcripts. They will be part of the cgranges query results when fragments are mapped to transcripts
- If after candidate pruning, the only remaining matches are nascent RNAs, then those nascent RNAs now have evidence of existence. Of course, genomic DNA also competes with nascent RNA, but that is a separate problem.



## Equivalence class sorting?

This seems to affect non-deterministic behavior of the EM. Is this necessary?

    // ---- Deterministic ordering ----
    // Multi-threaded BAM scanning produces fragments in non-deterministic
    // order.  The unordered_map above inherits that non-determinism in both
    // (a) the iteration order of equiv classes, and (b) the row order of
    // units within each class.  Since the EM E-step accumulates column sums
    // over rows, and FP addition is non-associative, different row orders
    // produce ULP-level differences that SQUAREM amplifies across iterations
    // potentially causing large cascading output differences.
    //
    // Fix: sort equiv classes by comp_idx, and sort rows within each class
    // by their log-likelihood fingerprint.  This makes the EM iteration
    // fully deterministic regardless of input fragment order.


## run_locus_em_native()

When/why is this function needed?


## C/C++ Optimization: AVX2

SIMD provides enormous speedup for the EM
We wrote a NEON specific fast_exp() implementation
We need an AVX2 implementation for linux

Eventually the production system will be Linux


## Edit distance

- when reads overlap, edit distance can get double counted in the overlap portion of the reads
- STAR avoids this with the 'nM' tag

- I think for most tools, we need to parse the MD tag of the BAM file, find the exact positions of the mismatches/edits, and then merge them between reads. It would require MD tag parsing and merging code. Sort of complicated.

- Prefer to just use STAR for now



## Mappability

### Genomic footprints

Genomic footprint is an important aspect of the coverage density model. It is the denominator upon which we compute coverage density. Without accounting for mappability, gDNA will be underestimated because the unmappable genomic regions will inflate the denominator.

Solution: as part of indexing, compute "mappability" tracks across the genome. 

Complicating factors:
- Mappability depends on read length
- Paired-end reads may "anchor" unmappable regions allowing them to be resolved.

Solutions:
- Compute mappability tracks for different read lengths
- Estimate total genomic footprint as a function of read length
- This will require mapping billions of reads. For every read length, we must map ~4 billion reads.


## Output 

- revised outputs (2026-03-22)
- need to review the output files carefully and see if changes are needed

## Post-EM Pruning

- pruning was recently removed
- now we have 'discrete' fragment count assignment as the default


## (RESOLVED) Fragment length distribution

- Review FL distribution observations
- Should be spliced fragments with unambiguous fragment length (sometimes fragment can map to different transcripts and have different fragment lengths)
- Need to ensure nascent RNAs are not interfering with FL distribution observations