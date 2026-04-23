# TODO


## Summary (output)

- Put "quantification" section towards the top with high-level results
- Fragment length distributions should be the last thing


## Calibration 

## Calibration: beta binomial versus binomial

- Current implementation is Binomial with strand floor and various other fixes
- Prior implementation used beta binomial which required estimating parameter kappa which was not straightforward and required limits


## Calibration: estimate_kappa_sym

This was computed at time where we were trying to constrain gDNA to be "symmetric" by strand so that gDNA strand assignments needed to approach 50/50 with some variance (binomial with overdispersion) and this was trying to assess the variance/overdispersion. We may not need this at all anymore.




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

