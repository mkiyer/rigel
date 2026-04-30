# TODO


## Calibration

Currently, calibration does a good job of obtaining high confidence gDNA fragments from which to train the gDNA FL distribution. However, the calibration system fails to accomplish its other main objective -- estimate locoregional gDNA levels within each genomic region. Locoregional gDNA estimates are what makes up the "prior" used in the EM. If we imagine a total amount of "prior" weight (c_base=10 for example), we wish to divide thus first between gDNA and total RNA (total RNA = mature RNA and nascent RNA). I think the first priority is to accurately estimate gDNA levels during calibration. We are making progress towards accurate gDNA estimation but we have not yet achieved this. Starting from the simple mini-genome single multi-exonic transcript case is helpful, but this is just the beginning. Currently, gDNA is estimated from intergenic, intronic, and exon-incompatible fragments. There is gDNA present as fragments completely contained within exons (exon-contained category) but these fragments are not utilized. Thus there is a population of gDNA that is not being estimated and is unaccounted for. Our strategy needs to estimate gDNA levels without being "blind" to certain regions of gDNA. Our intention is to eventually predict gDNA in all regions -- including the exon-contained regions where it is most difficult to predict. Having a reasonably accurate "pool level" estimate of gDNA is critical. First, we need to separate gDNA from RNA (RNA = mature RNA and nascent RNA). Then, we may need to distinguish nascent from mature RNA but that is a separate concern. So now we need to design a robust gDNA estimation system that will work for both unstranded and stranded libraries. Let's design it. Problem definition -- given a genomic region (e.g. "Locus") or set of genomic regions (e.g. a "MultiLocus") estimate the fraction of gDNA vs total RNA in the locus.

Fragments need to be categorized first:

Spliced fragments
- artifacts (splicing blacklist -> artifact category)
- mature RNA

Implicitly spliced fragments / non-contiguous alignment blocks
- intron within unsequenced gap between reads -> spliced equivalent
- other out-of-range fragments -> unannotated isoforms, chimeras, etc. -> held out of calibration

Unspliced fragments: aligned blocks to contiguous genomic region, no intron within unsequenced gap, fragment length within limit
- intergenic (gDNA only)
- overlapping intergenic AND transcript (gDNA vs extended transcript start/end)
- contained within transcript
-- entirely exon-contained (gdna + mature rna + nascent rna)
-- partially overlapping exon (gdna + nascent rna)
-- completely intronic (gdna + nascent RNA)


### Region "math"

Entire genome needs to be partitioned into contiguous non-overlapping "regions". 

Genomic region definition:
- ref, start, end
- tx_pos_bp, tx_neg_bp (bp of overlapping transcript on + and - strand)
- exon_pos_bp, exon_neg_bp (bp of overlapping exon on + and - strand)
- n_unspliced



## Implicit splicing

- Not handled yet


## Chimeric fragments

- How are these flowing through buffer?
- Phase 6



## Calibration

We are going to need to redesign the BAM scanning phase to correctly support our new calibration module

The main steps:

1) The first step towards estimating the gDNA FL is obtaining a reliable/robust pool of unspliced fragments.

2) The second step is to categorized unspliced fragments by exon and transcript overlap.

3) The third step is to consolidate fragments likely to be gDNA candidates and estimate the gDNA FL distribution


### BAM scanning phase should gather calibration data correctly

Each fragment has multiple aligned blocks. Each aligned block is a contiguous genomic region. The full genomic footprint of a fragment is the entire span of its aligned blocks.

We wish to ensure that the BAM scanning phase properly identifies "spliced" and "unspliced" fragments.

Spliced fragments may include:
- EXPLICIT: Fragments where reads have CIGAR "N" (skip) where the skip matches an annotated splice junction and passes the splice junction "blacklist" filter
- IMPLICIT: Fragments without an explicit CIGAR (skip) that have an unsequenced "gap" between aligned blocks (fragment length must be greater than 2 x read length, leaving an unsequenced gap in the fragment), where read1 and read2 match exons of the same transcript and the unsequenced "gap" contains an intron (splice) of that transcript. Presumably the splice junction was not explicitly sequenced but is present in the unsequenced gap.

Spliced fragments should NOT compete for genomic DNA.

Unspliced fragments include:
- TRUE: Fragments without CIGAR "N" (skip) AND (aligned blocks either overlap (no unsequenced gap) OR (aligned blocks span a contiguous genomic region without an "intron" in the unsequenced gap). The genomic span should be a single contiguous genomic region

- ARTIFACT: Fragment with CIGAR N Skip that is "blacklisted" because it is found in the splicing blacklist (the blacklist was created by aligning contiguous regions of DNA using a spliced aligner resulting in false positive spliced alignments). The aligned blocks from an artifact fragment should not be combined into a single genomic footprint. These fragments can be partitioned into distinct "blocks" for calibration purposes. One fragment -> multiple "blocks" for calibration where each block is an unspliced contiguous genomic region. It may be best to "hold out" these artifact fragments or consider them ambiguous for the first version of the calibration tool.

During the BAM scanning phase, fragments should be correctly classified and feed into the calibration step.

### Unspliced fragment pipeline

We consider unspliced, UNIQUELY ALIGNED fragments from the BAM scanning phase as initial candidates for gDNA.

Unspliced fragments should have defined, reliable contiguous genomic region span (ref, start, end) with a defined fragment length.

We then categorize unspliced fragments by their overlap with index features. Specifically overlap with transcripts and exons.

At this point each unspliced fragment has a defined fragment length!

The data need for categorization are:
- exon_overlap_pos: number of bases of fragment that overlaps exon from transcript on the positive strand
- exon_overlap_neg: number of bases of fragment that overlaps exon from transcript on the negative strand
- transcript_overlap_pos: number of bases of fragment that overlaps transcript on the positive strand
- transcript_overlap_neg: number of bases of fragment that overlaps transcript on the negative strand

Once we have this data, we can assign categories:

- INTERGENIC: ZERO overlap with transcripts (and by definition exons) on either strand
- INTRONIC: (ZERO overlap with exons on either strand) AND (non-zero overlap with introns).
- EXON_CONTAINED: ENTIRE fragment overlaps exons (fully contained within exons)
- EXON_INCOMPATIBLE: fragment partially but does not completely overlap exons.

We also assign each fragment a merged 'strand':
- none (undefined strand)
- pos (positive strand)
- neg (negative strand)
- ambiguous (overlaps transcripts on both strands)

The previous "region partition" github code (now removed) created a genome-wide interval index to allow rapid interval overlap query and efficient computation of exon/transcript overlap in base pairs. This is removed. We need to consider the most efficient solution to compute this in the current framework.


### Defined pool for gDNA FL estimation

The above steps (unspliced fragment pool, unspliced fragment categorization) are needed to find a pool of unspliced fragments that are compatible with gDNA.

We then build the initial gDNA FL model using the following categories:
1) INTERGENIC
2) INTRONIC
3) EXON_INCOMPATIBLE: here is where a "tolerance" in bp comes into play. If "overhanging" portion that does not overlap exon is within a tolerance (default 3bp) we will not include it in the gDNA FL model (at a small number of bp overhang, the fragment may simply be misaligned)

We then fit the gDNA FL model.

There is a 1D mixture EM step to refine the gDNA FL model using the RNA FL model and gDNA FL model.


## Summary (output)

- Put "quantification" section towards the top with high-level results
- Fragment length distributions should be the last thing




## cgranges interval searching

### possible to serialize cgranges index? saves time and memory...

### should we "lazily" build cgranges index if only used during one stage of pipeline? for example, calibration requires cgranges index for mappability but then might not be used again. could lazily build and free, potentially decreasing overall RSS




### Hybrid capture support

- If we are running hybrid capture, we need the calibration module to only use "captured" regions to train itself. We need to provide 'on target' and 'off target' regions.
- Defer for now



## Locus stats

Enable locus status to be written to the output


## Nascent RNA count set to zero? 

Copilot suggseted that synthetic nRNA counts are being set to zero? Is this some kind of gating?


### nascent RNA gating?

- formalize nascent RNA gating policy and implement this
- "gate" nascent RNA only when there are intronic reads. if only exonic reads, nascent RNA component not activated


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

