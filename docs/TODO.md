# TODO


## Fragment length distribution

- Review FL distribution observations
- Should be spliced fragments with unambiguous fragment length (sometimes fragment can map to different transcripts and have different fragment lengths)
- Need to ensure nascent RNAs are not interfering with FL distribution observations


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
