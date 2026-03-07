# TODO




## Benchmark

- benchmark.py
- salmon needs a salmon index (full transcriptome)
- kallisto needs a kallisto index (full transcriptome)
- build or point to it


## Fragment length distribution integration

- this is planned
- strategy to weight observations by uncertainty
- build RNA and gDNA fragment length distributions
- use distributions to improve deconvolution


## Edit distance

- when reads overlap, edit distance can get double counted in the overlap portion of the reads
- STAR avoids this with the 'nM' tag

- I think for most tools, we need to parse the MD tag of the BAM file, find the exact positions of the mismatches/edits, and then merge them between reads. It would require MD tag parsing and merging code. Sort of complicated.

- Prefer to just use STAR for now



## Overhang based gating

1) We should review our overhang-based GATING (winner take all) -- what information is needed there? Do we need our strand and fragment length likelihood for overhang gating? It seems like we do.. because if there are fragment-transcript matches on both strands (overlapping antisense transcripts) the fragment could overhang by different amounts but the strand likelihood REALLY matters for gating.. we wouldn't want to gate solely based on overhang alone because we would potentially throw out the more likely matches.. this only happens in when there are overlapping opposite direction transcripts though.. if same strand transcripts we can probably gate more easily... that leads me to an idea -- For a SUBSET of fragments (they need to be unambiguous with respect to strand.. e.g. isoforms on same strand) we might be able to do overhang gating earlier.. at time of fragment resolution.. but again.. i am starting to worry about our gating strategy because we DO want strand, fragment lenght, and overhang to all play into this. it's possible that transcripts with some overhang could score better if you factor in all of the likelihoods (strand, frag len, overhang, etc) and so gating on overhang may be dangerous. 


## Fragment length model

- gDNA and RNA should have different models
- different models should be used to help predict
- output shouldn't necessarily report exonic, intronic, intergenic, etc. that's just for initialization


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



### False positive splice junctions

We can "blacklist" splice junction artifacts by doing genome-wide mapping. Align to the genome and detect cases where gDNA sequence -> spliced alignment (artifact). Compile a "blacklist" from these cases.




## Output 

- rename 'tpm' -> 'mrna_tpm' and move next to mrna columns
- for tsv files, truncate floating point for readability






## 2026-03-03 Benchmarking

Equivalence class collapse (High severity) — 402/1522 transcripts (26%) across 9/10 regions receive identical estimates because they share the same fragment compatibility pattern. The EM splits counts evenly. Salmon resolves some of these via positional bias modeling.

Gene-level phantom smearing (Medium) — The EM leaks sub-integer fractional counts to zero-truth genes (e.g., 0.01 counts), causing rank shifts of up to 10 positions on Spearman despite negligible absolute errors. Fix: post-EM thresholding at ~0.5 counts.


Three transcript-level losses (Low) — MALAT1_NEAT1, TTN, EGFR all trace back to equivalence class collapse (29-42% collapse rates).


## Post-EM Pruning

Currently we use evidence ratio (percent count data relative to the prior) with recommendation of 10%.

Another option would be just zeroing out transcripts with < 0.5 counts and running one final EM step. Should we do the zeroing by absolute counts or by evidence ratio?


## Index

Index takes a long time to load?




## (RESOLVED) 2026-03-03 nRNA leak bug

False nascent RNA classification (Low) — Up to 33/100k reads misclassified as nascent RNA in pristine data, mainly at HOXA_cluster (short introns create ambiguity).

nRNA leak increased at HOXA_cluster (33→108 reads) — this is likely due to different abundance draws rather than a pruning effect.

It appears that a very small amount of nRNA is leaking when it shouldn't. This is happening in "prisine" scenarios with nRNA set to 0.0 and "oracle" reads that are perfect. If no alignment errors and no nRNA, any nRNA we observe is likely a bug. We now have a locus that we can analyze (HOXA) and the annotate BAM file strategy that will give us insight into which fragments are being inappropriate misclassified as nRNA.

- This was due to a global prior of 0.5 applied to the 'eta' variable
- Now we estimate the global prior from the data


## (RESOLVED) Nascent RNA linked to mRNA

2026-03-03: Created a new model / parameterization based on the kinetic model of transcription -> nRNA -> splicing -> mRNA -> degradation -> null. This model intrinsically links nRNA with mRNA. Previously, mRNA and nRNA were solved for independently by the EM.

## (RESOLVED) Unstranded RNA-Seq

2026-03-03: For unstranded data, we cannot rely on strand specificity to estimate gDNA vs RNA. The reads are evenly distributed. In this case, we can still use coverage density to estimate genomic DNA.
- intergenic: gDNA
- intronic: nRNA + gDNA
- exonic: mRNA + nRNA + gRNA
- spliced (annotated): mRNA
- spliced (not annotated): mRNA or mapping artifact

The current solution is to create a hybrid estimator that uses coverage density and strand specificity to estimate the fraction of gDNA vs nRNA vs mRNA and initialize the EM.

This uses:
W_strand = (2*ss - 1)^2 
W_density = 1 - W_strand

The strand-specificity is the "knob" controlling whether density or strand information is used. When ss = 0.5, W_strand = 0 and W_density = 1.


## (RESOLVED) Empirical bayes shrinkage of priors for "eta" = fraction of nascent RNA / total RNA

Same 4-level hierarchy (transcript → TSS → locus-strand → weak)

- This is implemented as a "cascading" if else statment
- Can this be done with a "smooth" approach that incorporates all sources of information in the hierarchy? The if-else statement seems "choppy" and wrong.

2026-03-03: Shifted to empirical bayes. 


## (RESOLVED) Unambiguous overlap in HulkIndex

Added the intervals of 'unambiguous introns' into the index to make overlap queries more efficient


## (RESOLVED) BAM parsing XS tag bug

I want to ensure robust BAM tag parsing. Do you need to look at examples and learn best practices?

Read htslib spec here:
https://github.com/samtools/htslib/blob/develop/htslib/sam.h

Learn how pysam parses tags? Or look at samtools code?
https://github.com/pysam-developers/pysam/tree/master/pysam

If not done already, please add unit tests to ensure robust tag parsing Do we need a C/C++ unit testing framework to directly support testing the C/C++ implementation? Google Test?

Here is a the strand parsing logic. When a read is spliced, its BAM CIGAR string will contain an 'N' (skip). We can obtain the precise splice junction coordinates (the intron) by parsing the read and the CIGAR string. At this point, we 'could' look for annotated splice junctions in our HulkIndex by matched (ref, start, end, <strand>) but would have to query the index for both POS and NEG strands. But even without an XS tag we could lookup ANNOTATED splice junctions by checking our index and both strands. However, the aligners are supposed to detect unannotated splice junctions as well and will use the DNA splice motifs (GT/AG, etc) to mark the strand as POS or NEG. For STAR, the XS tag is the strand of the matching splice junction. However, for minimap2, this is a flag that tells us whether we need to flip or not flip our read. So our code should be able to parse splice junction with an XS tag (STAR), a 'ts' tag (minimap2) requiring different logic, or NO tag (in which case we could query our annotated splice junction index and obtain the strand from matches). Once we have the strand of the splice junction, we can store it in a list of splice junctions for that read. Then we MERGE reads from multiple BAM records. (for pe reads, supplementary reads, etc). Merging the splice junctions should be more or less a set union procedure, as the splice junctions are exact. If reads have splice junctions on opposite strands (read 1 + and read 2 -) we have a major problem -- it indicates an ambiguous read which is really only possible if there is a chimera (genome alteration or trans-splicing) creating unnatural/abnormal transcripts. We are SIMULATING perfect reads on a perfect genome and should never read ambiguous strands (read1 and read2 have splice junctions on opposite strands). This would be a bug.


## (RESOLVED) Refactor 'eta' → 'nrna_frac'

Renamed all occurrences of 'eta' / η to 'nrna_frac' across the entire codebase:
Python source, C++ em_solver, CLI flags (--nrna-frac-*), tests, and docs.
708/708 tests passing after rename.

nrna_frac = nascent RNA fraction = (nRNA / (total RNA))

Hierarchy: nrna_frac_t, nrna_frac_tss, nrna_frac_locus, nrna_frac_global



## (RESOLVED) Refactor 'unique' → 'unambig'

Renamed EM-resolution "unique" identifiers to "unambig" throughout the codebase
to distinguish from "unique mapper" (alignment context, NH=1).

Convention:
- `unambig` / `ambig` — fragment-transcript resolution (one candidate hit vs. multiple)
- `unique` / `multimapping` — read alignment (one genomic location vs. multiple)

Renamed: FRAG_UNIQUE→FRAG_UNAMBIG, unique_counts→unambig_counts,
unique_totals→unambig_totals, assign_unique→assign_unambig,
deterministic_unique_units→deterministic_unambig_units,
em_routed_unique_units→em_routed_unambig_units, n_gdna_unique→n_gdna_unambig,
frag class label "unique"→"unambig", C++ unique_totals→unambig_totals.
708/708 tests passing after rename.


## (RESOLVED) Strand Model

- StrandModel "hard coded" cutoff of # of reads
- output shouldn't necessarily report exonic, intronic, intergenic, etc. that's just for initialization
- really just two models, RNA and gDNA

Smoothing towards 0.5 with pseudocount parameter
kappa = 2 default

2026-03-03: only use spliced reads to train. other reads are contaminated.


## (RESOLVED) Incorrect nRNA and gDNA initialization

Change to a 'linked' total RNA = nRNA + mRNA model led to the creation of a subtle initialization bug where gDNA and nRNA reads were being double counted. This has been resolved.


## (RESOLVED) Any vestigial / stale hard cutoffs for strand specificity?

Previous versions of the code had hard cutoffs for strand specificity. Do those still exist? These should be removed because our new initialization scheme uses both coverage and strand specificity together. We should not longer require cutoffs for strand specificity anywhere in the code. Do a deep dive thorough search for any remaining strand specificity conditional statements based on some cutoff. 

## (RESOLVED) pthread vs openmp

can we simplify by using pthread for our parallel locus EM? what is the benefit to openmp?

2026-03-06: Replaced OpenMP with std::thread. The locus EM now uses a simple atomic-counter work-stealing loop (chunks of 16 loci) with std::thread workers and std::atomic<double> for the single shared accumulator. Removed libomp dependency from CMakeLists.txt and mamba_env.yaml. All 766 tests pass.
