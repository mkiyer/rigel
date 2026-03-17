# TODO


## Reference ID (not string)

Reference names are stored throughout the code base as strings. We should change this to use their integer IDs for efficiency.


## gDNA siphoning problem

Currently gDNA siphons mRNA + nRNA fragments. Initialization is a major problem. The EM solver might also need constraints (strand "symmetry" constraints where gDNA symmetry is 50% sense 50% antisense and RNA symmetry is SS / (1-SS).


## C/C++ Optimization: AVX2

SIMD provides enormous speedup for the EM
We wrote a NEON specific fast_exp() implementation
We need an AVX2 implementation for linux

Eventually the production system will be Linux


## Equivalence class sorting?

This seems to affect non-deterministic behavior of the EM. Is this necessary?




## empirical bayes

- review gdna empirical bayes setup which requires a bunch of constants

- review nascent RNA 'eta' empirical bayes setup


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


