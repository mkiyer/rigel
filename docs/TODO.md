# TODO


## gdna capture weighting

I'm just thinking out loud here.

We already have regions partitioned.
- Intergenic (gdna)
- Intronic (gdna + nrna)
- Exon-intron (gdna + nrna)
- Exon (gdna + nrna + mrna)

Given that nascent RNA tends to be very sparse genome-wide, intergenic + intronic + exon-intron are reasonable to approximate gDNA density.

Each region is a genomic interval and after this step we have gDNA density estimates for each region. EXON regions are also estimated by projection of the bordering exon-intron regions. So all regions should have gDNA estimates. That should already be happening in the code.

The simplest next step is to figure out how to weight each region. Could each region could be assigned a weight?

What if we divide/partition the total gDNA "mass" among the regions? Each region gets a weight proportionate to its contribution to the total gDNA mass.

So for example let's say there is the 100kb locus with 1kb exon and 99kb intron. We have the exon gdna densities (computed from the exon-intron boundaries) and intron gdna densities. There are 1000  total gDNA fragments in the 100kb locus, but 950 of them are consolidated on the exons (950/1kb = 0.95 frags/bp density) and 50 of them are in the introns (50/99kb = 0.0005555 frags/bp density).

We then would assign per-region weights such that the GLOBAL sum of the weights equals 1.0, and the regions are then weighted by their contribution to the total global gDNA mass.

This seems extremely simple to me. The aggregate gDNA locus level estimates are then multiplied by their weights so that the high-density regions contribute much more weight than the low-density regions.

I do recognize that there is a "circular" logic issue here where we are weighting gDNA estimates by their densities.. kind of but not exactly effectively using density squared as the weight rather than density alone. On that note, raising the density to a power effectively squashes small values and retains larger values, effectively distributing density across more highly-weighted regions.

You mentioned an idea to subsample the total region population to build a density distribution so that we are using a portion of the total regions to construct a density distribution and then using that distribution predict per-region weights. I'm not sure if that helps but repeated subsampling and computing weights would incorporate some aspect of uncertainty in the region measurements.

Speaking of uncertainty, we already had a plan to enhance how we compute bayesian gdna prior. The current bayesian gdna prior does not include uncertainty in any way. It is computed as a flat pseudocount that adds. The gDNA estimation has an uncertainty component that ideally should be modeled.

Uncertainty comes from several sources. First, tiny regions have much greater uncertainty than large regions (big difference in exposure). Exon-intron boundary regions are small and there will be substantial uncertainty there.

Uncertainty comes from not-modeled factors that create gDNA variation such as a GC content.

When we apply gDNA weighting, we should consider incorporating uncertainty modeling at the same time.







We don't have enough information yet to be able to change our code or parameters. We need to delve deeper into this issue.  FIRST, I think we can isolate the problem by eliminating/discarding multi-mapping reads. First audit the code, there should be a CLI parameter include multimapping that turns on/off inclusion of multi-mapping reads. Interrogate the code flow and ensure that this is implemented correctly. IF we turn OFF multimapping reads, we should break the mega-locus issue. If the mega-locus problem is the main contributor to gDNA -> RNA errors, we should be able to diagnose and visualize this after discarding multimapping reads. Without multimapping reads, we should see a huge performance boost. Please do this first. 

## ambiguous region errors

We don't have enough information yet to be able to change our code or parameters. We need to delve deeper into this issue. I want you to find regions where gDNA fragments are mispredicted to be RNA using the annotated BAM output file. Select up to 5 regions that are particularly egregious. Show that these regions are so we can look at the BAM file on the genome browser and a get a sense of what the reads look like -- find precise locations where there are a very high number of mispredicted fragments. Analyze these on an individual per-fragment level to understand what is happening to the likelihoods of the fragments. We need to find out exactly what is happening to these fragments. 


## mappability corrected effective length
   I agree strongly with your concern. Intergenic calibration without mappability-corrected effective length is fragile, and hybrid capture breaks naive intergenic extrapolation. The denominator should become something like accessible/capturable effective length, not raw genomic span. For capture libraries, gDNA exposure should include probe/bait targetability or an empirically learned capture profile.

## gdna prior uncertainty
  **Treat the gDNA prior as uncertain, not as a hard pseudocount.**
   Right now the calibration-derived `gdna_prior_count` enters as a physical Dirichlet count. In unstranded data, especially when calibration channels may contain nascent RNA, that prior should be downweighted by an evidence-quality factor. Good knobs to test: cap per-locus `eta_g / n_em_fragments`, shrink the prior toward zero when strand contrast is unidentifiable, and propagate calibration uncertainty instead of using a point estimate.

## gdna vs nrna eff length normalization
   Current gDNA length normalization is not totally unconstrained; it is baked into the gDNA log likelihood using a local sampling window based on transcript span plus a gDNA flank. Synthetic nRNA then gets transcript-like effective-length normalization in EM. For long gene spans, the nRNA “shorter than gDNA” advantage is probably tiny, and for exon-contained fragments mRNA has the real effective-length advantage. I would audit whether gDNA should use the full locus/capture-accessible union instead of an anchor transcript span, but avoid arbitrary penalties that create false nRNA.


## Remove 'coverage_weights'?

Not used as prior anymore
Used as warm start
But does it do anything?


## Rigel index with GTF duplicates

GENCODE GTF has duplicates (why?)
Need option to drop duplicates during index build
Keep lexicographically smallest transcript id?

## Performance optimization


## Locoregional calibration

Currently, calibration collect global information to estimate the gDNA prior for bayesian EM.  This might be improved by incorporating locoregional information, which may reflect copy number changes.

Given a particular Locus (genomic region containing transcripts) this approach would compute gDNA densities from the flanking regions of some genomic distance (how much? 10kb, 100kb? 1Mb?). The locoregional gDNA estimates would shrink to the global estimates using empirical bayes shrinkage and then feed into the bayesian prior for the Locus of interest.



# Polars migration

Polars is apparently much faster than pandas. Might be ideal for rigel

## cgranges interval searching

### possible to serialize cgranges index? saves time and memory...

### should we "lazily" build cgranges index if only used during one stage of pipeline? for example, calibration requires cgranges index for mappability but then might not be used again. could lazily build and free, potentially decreasing overall RSS







## Gene counting

Ensure gene quantification is correct.

Regarding gene quantification. Yes, we can sum counts of the transcript isoforms of the gene. Annotated single-exon transcripts are associated with a gene and will be included in those counts (they are both mature and nascent so this is okay). For gene quantification, we MUST exclude synthetic nascent RNA (synthetic nascent RNA should not be associated with a gene. Remember, we can sum counts, but when we compute transcripts-per-milion (TPM), we use a weighted average of counts respecting different effective lengths.




## Transcript versus pseudogene multimapping (edit distance improvement)

It appears that transcript vs processed pseudogene identifiability is the largest remaining source of error in rigel. When running with minimap2, rigel fails to distinguish the true transcript from its processed pseudogene, resulting in massive counting error.

Our hope is that the 'true' transcript will contain fewer mismatches/edits and can be identified by superior alignment quality compared to other competing alignments.

We wish to improve how we distinguish multimapping alignments. We want to be able to prioritize the 'true' alignment relative to 'false' alignments.

The rigel tool currently tries to do this by parsing the 'NM' tag in the BAM file. This tag is supposed to represent the 'edit distance' between the read and the reference. Using 'NM' is intended to help distinguish alignments by the number of edits between the read and the reference. However, some aligners may or may not be setting the 'NM' tag appropriately. Some tools only set 'NM' when there are insertions/deletions. It is unclear how minimap2 sets the 'NM' tag.

The 'MD' tag in SAM/BAM files is a companion to the CIGAR string and is intended to identify the exact positions of the mismatches/edits. 

Here's an excerpt from the PDF describing SAM/BAM optional fields, that defines and explains the MD tag:

### MD tag definition

MD:Z:[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*

String encoding mismatched and deleted reference bases, used in conjunction with the CIGAR and SEQ fields to reconstruct the bases of the reference sequence interval to which the alignment has been
mapped. This can enable variant calling without requiring access to the entire original reference. The MD string consists of the following items, concatenated without additional delimiter characters:
- [0-9]+, indicating a run of reference bases that are identical to the corresponding SEQ bases;
- [A-Z], identifying a single reference base that differs from the SEQ base aligned at that position;
- \^[A-Z]+, identifying a run of reference bases that have been deleted in the alignment.

As shown in the complete regular expression above, numbers alternate with the other items. Thus if two mismatches or deletions are adjacent without a run of identical bases between them, a ‘0’ (indicating a 0-length run) must be used to separate them in the MD string. Clipping, padding, reference skips, and insertions (‘H’, ‘S’, ‘P’, ‘N’, and ‘I’ CIGAR operations) are not represented in the MD string. When reconstructing the reference sequence, inserted and soft-clipped SEQ bases are omitted as determined by tracking ‘I’ and ‘S’ operations in the CIGAR string. (If the CIGAR string contains ‘N’ operations, then the corresponding skipped parts of the reference sequence cannot be reconstructed.)
For example, a string ‘10A5^AC6’ means from the leftmost reference base in the alignment, there are 10 matches followed by an A on the reference which is different from the aligned read base; the next 5
reference bases are matches followed by a 2bp deletion from the reference; the deleted sequence is AC; the last 6 bases are matches.

### Double-counting edit distance when paired-end reads overlap

- When a fragment size is less than 2X read length, then the paired-end reads sequence a common portion of the fragment. For example if the fragment length is 250bp and we do 2 x 150 bp paired-end sequencing, then the two reads redundantly sequence a common 50bp region of the fragment. When we calculate edit distance, we have to ensure that mismatches/edits in this common "double sequenced" region do not get double counted.

- The STAR aligner avoids this double counting with the 'nM' tag
- Minimap2 and other aligners do not appear to handle this issue


### Proposed solution and implementation steps

- Instead of using the annotated BAM file, let's start by creating a simulated scenario using the GAPDH gene. 
- We have a wonderful 'sim.py' script in 'scripts'. We can set this up to produce reads to GAPDH transcripts (and zero transcripts everywhere else). Then we can generate a tiny simulated scenario that can be used for debugging.
- The simulation needs the full human genome reference so that reads multimap to other locations.
- Generate simulated paired-end reads to transcripts from the GAPDH gene. GAPDH has many pseudogenes and reads often multimap. The current rigel tool with minimap2 has massive errors related to GAPDH abundance estimation.
- Create a simulated "oracle" BAM file to GAPDH
- Align the reads using minimap2
- Run rigel on the oracle BAM 
- Run rigel on the minimap2 BAM
- Trace the reads/fragments through both the oracle BAM and the minimap2 BAM.
- Find fragments that multimap. Understand how these are being processed. Understand how the likelihoods are calculated. Understand at a deep level why Rigel cannot correctly distinguish these reads.

#### Big question: can we distinguish transcript from pseudogene using any aspect of the alignments?

- NM tag?
- MD tag?
- Other alignment results?
- How would we design this?



## Genome to transcript coordinate search

This code may be duplicated/redundant in more than one place? It is being used in scoring and also in fragment length calculation. Can we consolidate code? Can we improve the efficiency of data structures?


## Unannotated splice junctions

A fraction of splice junctions are "unannotated" in that they don't have exact matches to the reference. These can be either artifacts (aligner errors) or 'novel' isoforms/transcripts. There is a lot of information that we are not taking advantage of yet:

- partially annotated: 3' site, 5' site, or both
- (easy) is the 5' splice site annotated
- (easy) is the 3' splice site annotated
- (medium) how much 'anchor' on either side of the intron?
- (hard) does this junction match any of the 'blacklist' splice junctions -- requires building a splice junction blacklist.






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

