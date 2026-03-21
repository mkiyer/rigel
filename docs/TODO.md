# TODO


## Post-EM discrete abundance estimation

As the EM is running we maintain posterior probabilities and fractional transcript counts. However, the input data are discrete molecular sequences, where each molecule can only arise from *one* transcript/nascent RNA/gDNA in the genome even there may be more than one possible source. If we are to stay true to the data, then our final output should be in the form of discrete counts rather than fractional counts.

We *definitely* should not discretize the EM algorithm itself.

My idea is that after the EM converges, we should do an additional iteration over the fragments. Post-EM, each fragment has a set of compatible transcript/nRNA/gDNA candidates each with a posterior probability. We now *know* the probabilities that the fragment could have arisen from each of its possible compatible hits. So, what I propose we is as we iterate through the fragments again post-EM, we can assign each fragment to a *single* candidate transcript/nRNA/gDNA. There are a few different policies we could use for assignment: 1) always assign to the highest probability candidate (e.g. "greedy" or "best", or "winner take all" policy), 2) "multinomial" mode - assign fragment with multinomial sampling based on the posterior probabilities (ex. fragment compatible with transcript A p=0.05, transcript B p=0.8, or transcript C p=0.15, we then 'randomly' assign the fragment to *one* of its compatible transcripts based on the probabilities.. 85% chance it gets assigned to A, 15% it gets assigned to C, and 5% change it gets assigned to B). 

If we assign fragments using 'multinomial' sampling, we can also consider 'pruning' from the list of compatible transcripts any transcript with posterior probability that is less than a threshold.. so that very very low probability candidates are eliminated from being selected altogether.

Please evaluate this idea. Do you think it is a good idea? Is it possible to implement? What would the implementation look like? What would be required? This would introduce parameters: --post-em-sampling ('best', 'multinomial'), --post-em-prob-threshold (default 0.05). You might suggest different names that work better or are more appropriate. Please give me your ideas and a plan for moving forward (no implementation yet).


## "Annotated" nascent RNAs

Reference transcript annotations vary considerably throughout the human genome. I found that there are a subset of transcripts that effectively already have nascent RNAs annotated.

For example, the MALAT1 long non-coding RNA has many spliced isoforms, but it also has a single-exon transcript isoform that spans the entire gene and "contains" all of the spliced isoforms (it overlaps the spliced transcripts and 'contains' them). In this example, the 'nascent RNA' already exists in the reference annotation!

An assumption of nascent RNA modeling is that nascent RNAs do not already exist as transcripts in the reference gene annotation. When the 'nascent' transcript is *annotated* in the reference, nascent RNA modeling will be affected. First, single-exon genes by definition cannot be modeled as nascent RNAs because there is no way to differentiate mature from nascent single. Second, adding a nascent RNA 'shadow' transcript to a gene/locus that already effectively has annotated nascent RNA just adds redundancy.

To solve this problem we need a few things. 

- First, we need a way to flag annotated reference transcripts as 'nascent RNA equivalent' transcripts, so that they can potentially be treated as such.
- Second, we need to merge or consolidate nascent RNAs with reference annotated nRNA equivalents so that they coalesce intead of competing. 

Properties of nascent RNA:
- single exon
- overlaps an entire transcript span from start to end

An extension of this issue:
- Currently, every combination of transcript start/stop gets its own nascent RNA shadow. Sometimes, transcript starts can vary by only a few bases. Since the transcripts aren't 'exactly' the same start/stop, they will emit distinct nascent RNA shadows, even thought the nascent RNAs may be almost completely identical (except for a few bases at the start or end).


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



## dead/scale code

— Unspliced sense/antisense: Kept as-is (C++ fused_score_buffer requires the arrays in its signature; removing would be too invasive for negligible benefit)

- intronic_sense/antisense accumulators are populated by the C++ scoring pass but are not used. These are likely dead/stale code.



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


