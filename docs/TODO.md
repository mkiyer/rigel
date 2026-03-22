# TODO


## Nascent RNA gating?

- We could require evidence and gate nascent RNAs
- Evidence would be fragment within an 'unambiguously intronic' genomic interval
- In the new model, nascent RNAs are just (synthetic) transcripts. They will be part of the cgranges query results when fragments are mapped to transcripts
- If after 'winner take all' fragment resolution, the only remaining matches are nascent RNAs, then those nascent RNAs now have evidence of existence. Of course, genomic DNA also competes with nascent RNA, but that is a separate problem.

## Winner take all fragment resolution

- This is a mechanism to reduce the number of compatible transcripts for a fragment.
- A fragment may be "compatible" with many transcripts, but some of the transcripts are much better candidates than others. Instead of polluting the equivalence class (the set of compatible transcripts) with a bunch of terrible candidates, we have a step that prunes out the worst candidates and keeps the best.
- We need to investigate the mechanism for this equivalence class gating step. The original implementation is based on 'overhang' -- Overhang past a splice junction into the intronic space is heavily penalized, and overhang that extends outside the transcript boundary (transcript start site and transcript end site) is heavily penalized. The "winner take all" approach finds the "best" transcripts (transcripts with the least amount of overhang) from the list of candidates.

## Equivalence class sorting


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


## Equivalence class sorting?

This seems to affect non-deterministic behavior of the EM. Is this necessary?





## Edit distance

- when reads overlap, edit distance can get double counted in the overlap portion of the reads
- STAR avoids this with the 'nM' tag

- I think for most tools, we need to parse the MD tag of the BAM file, find the exact positions of the mismatches/edits, and then merge them between reads. It would require MD tag parsing and merging code. Sort of complicated.

- Prefer to just use STAR for now



## Overhang based gating

1) We should review our overhang-based GATING (winner take all) -- what information is needed there? Do we need our strand and fragment length likelihood for overhang gating? It seems like we do.. because if there are fragment-transcript matches on both strands (overlapping antisense transcripts) the fragment could overhang by different amounts but the strand likelihood REALLY matters for gating.. we wouldn't want to gate solely based on overhang alone because we would potentially throw out the more likely matches.. this only happens in when there are overlapping opposite direction transcripts though.. if same strand transcripts we can probably gate more easily... that leads me to an idea -- For a SUBSET of fragments (they need to be unambiguous with respect to strand.. e.g. isoforms on same strand) we might be able to do overhang gating earlier.. at time of fragment resolution.. but again.. i am starting to worry about our gating strategy because we DO want strand, fragment lenght, and overhang to all play into this. it's possible that transcripts with some overhang could score better if you factor in all of the likelihoods (strand, frag len, overhang, etc) and so gating on overhang may be dangerous. 



## Unannotated splice junctions

A fraction of splice junctions are "unannotated" in that they don't have exact matches to the reference. These can be either artifacts (aligner errors) or 'novel' isoforms/transcripts. There is a lot of information that we are not taking advantage of yet:

- (easy) is the 5' splice site annotated
- (easy) is the 3' splice site annotated
- (medium) how much 'anchor' on either side of the intron?
- (hard) does this junction match any of the 'blacklist' splice junctions -- requires building a splice junction blacklist.


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


### Unannotated splice junctions

- partially annotated: 3' site, 5' site, or both


### False positive splice junctions

We can "blacklist" splice junction artifacts by doing genome-wide mapping. Align to the genome and detect cases where gDNA sequence -> spliced alignment (artifact). Compile a "blacklist" from these cases.


## Output 

- we still need to revamp the tool output
- current tool output does not include many of the recent updates or fixes (nascent RNA, etc)
- rename 'tpm' -> 'mrna_tpm' and move next to mrna columns
- for tsv files, truncate floating point for readability
- nascent rna output


## Post-EM Pruning

Currently we use evidence ratio (percent count data relative to the prior) with recommendation of 10%.

Another option would be just zeroing out transcripts with < 0.5 counts and running one final EM step. Should we do the zeroing by absolute counts or by evidence ratio?


