# TODO



## Calibration 


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




## Overhang likelihood penalty

It appears to be fairly common the minimap2 alignments that *should* span introns do not because they would have a 'tiny' anchor on one of the two exons. It seems that the aligner aligns the reads as unspliced with an intronic overhang instead of spliced with a CIGAR 'N' operation and a tiny anchor on the following exon. This is fundamentally an aligner problem.

The question here is whether we can set the overhang penalty to a value that penalizes but does not destroy the true candidate from among the different candidates.

The solution might *seem* easy when there is zero nascent RNA or genomic DNA contamination. The problem is that when genomic DNA contmination and nascent RNA are high, we have a major problem deconvoluting unspliced fragments (the fragment could be mature RNA, nascent RNA, or genomic DNA). Furthermore, the exonic overlap does not help, because we often have hybrid capture libraries that capture exonic nucleic acid and we enrich for DNA that way.

My gut feeling tells me that overhang penalty is too strict. However, the overhang penalty is fundamentally related to the level of genomic DNA contamination and nascent RNA. When genomic DNA and nascent RNA levels are low, then a small overhang should not incur a massive penalty. A penalty of 0.01 is very steep and hard to overcome. However, when genomic DNA levels are very high, we need to be more strict with overhang penalty, because a small amount of overhang may be a major signal discriminating gDNA or nRNA from spliced (mature) RNA.

Now that we have our gDNA calibration results, we are estimating gDNA levels in the sample. I wonder if we can set the overhang penalty to something more reasonable, possibly related to gDNA levels. However, we don't have nascent RNA levels as part of our calibration. Nascent RNA cannot be estimated globally as nascent RNA varies widely from transcript to transcript.

- Try running the tool with different settings for the overhang penalty. This should confirm the suspicion that tiny overhangs effectively 'kill' the correct transcript.

### Investigation results (2026-03-23)

**Status:** Wiring is already complete and correct. The penalty is configurable via:
- CLI: `--overhang-alpha` (0–1 linear penalty, converted to log internally)
- YAML config: `scoring.overhang_log_penalty`
- Default: `alpha=0.01` → `log(0.01) = -4.6052` per overhang base

**GAPDH case study using minimap2 BAM (50k simulated fragments, ss=0.95, gdna=0):**

The 1,656 "dominated" units (3.6% of all GAPDH EM units) were confirmed to be caused by **exactly 1 overhang base** in MANE vs 0 in the competitor (`ENST00000466525.1`).  The raw buffer shows `exon_bp(MANE) = read_length − 1` while the competitor has `exon_bp = read_length`. Source: `scripts/debug/confirm_oh_hypothesis.py`.

Overhang alpha sweep results (MANE ratio = mm2_TPM / oracle_TPM):

| alpha | oh_pen/base | MANE ratio | ENST396859 ratio |
|-------|-------------|------------|-----------------|
| 0.01  | −4.61       | 0.589      | 6.18×           |
| 0.10  | −2.30       | 0.597      | 6.52×           |
| 0.30  | −1.20       | **0.607**  | 6.71×           |
| 0.50  | −0.69       | 0.603      | 6.67×           |
| 1.00  | 0.00        | 0.471      | 10.63×          |

**Key conclusions:**
1. The overhang penalty has only marginal effect on MANE accuracy (+0.02 ratio across the full sweep).
2. **The dominant error is not the overhang penalty** — it is the **96.3% perfectly tied EM units** caused by shared exon architecture across all 11 GAPDH isoforms. The EM has no likelihood signal to distinguish MANE from other isoforms for these tied fragments.
3. Removing the penalty entirely (alpha=1.0) makes things **worse**: ENST396859.5 explodes to 10.6× over-count and MANE drops to 0.471. The overhang penalty actually provides a small but real filtering signal.
4. The adaptive approach (linking overhang penalty to calibrated gDNA fraction) is still worth exploring, but will not solve the fundamental isoform identifiability problem for within-gene transcript confusion.

**Remaining open question:** Why do ENST396859.5 (6×) and ENST619601.1 (5×) absorb so many MANE fragments? Length alone cannot explain it (ENST396859 is only 29 bp shorter). The likely cause is pseudogene multimapper fragments that land on exons shared between MANE and these shorter isoforms but not on MANE's unique 5′ region.

### Root cause — confirmed (2026-03-23)

Using oracle vs minimap2 EM unit classification (`scripts/debug/compare_em_classifications.py`, `scripts/debug/trace_alignment_mechanism.py`, `scripts/debug/junction_mismatch_analysis.py`):

**Oracle also has 98.9% tied units — but oracle quantifies correctly.** The tied units hypothesis was WRONG as the explanation for the oracle/mm2 divergence.

The real difference is:

| Category | Oracle | Minimap2 | Δ |
|----------|--------|----------|---|
| mane_absent (MANE not a candidate) | 260 | 2,919 | **+2,659** |
| mane_dominated (MANE loses to competitor) | 6 | 1,656 | **+1,650** |

These +4,309 "wrong seed" fragments bootstrap the EM to give higher weight to shorter isoforms. EM positive feedback then amplifies this seeding error into 6× over-count for ENST396859.5.

#### The physical mechanism: short exon anchor failure

MANE has a **52bp exon** at chr12:6534809–6534861, flanked by:
- a 240bp intron upstream: donor=6534569, acceptor=6534809
- the large 1632bp intron downstream: donor=6534861, acceptor=6536493

In **oracle BAM**: reads correctly produce 3-intron spliced alignments: `...M 240N 52M 1632N M...`

In **minimap2 BAM**: reads at this exon junction **fail to anchor the short 52bp exon**. The 240bp intron is skipped, producing: `13S 55M 1632N ...M` — the 13 soft-clipped bases are the "previous exon" anchor that minimap2 couldn't place, and the R2 exon-1 block of 55bp maps 3bp INTO the intron before the small exon.

Consequence for scoring:
- **MANE**: its exon starts at 6534809, so the mapped block [6534806–6534861] has only 52/55bp within MANE's exon → **oh=3** → penalty = 3×log(0.01) = −13.8 log units → MANE's posterior collapses → **MANE pruned from candidates**
- **ENST396859.5 / ENST466525.1**: these isoforms do NOT have the 240bp intron — their exon extends continuously through the [6534569–6534861] region → the full 55bp of the R2 block is within their exon → **oh=0** → they win

This is confirmed by junction comparison: **60% of mane-source mane_absent reads have DIFFERENT junctions in mm2 vs oracle**, always missing the (6534569, 6534809, 240) junction. And 99.2% of mane_absent units have MANE completely absent from the cgranges buffer (not pruned — genuinely not found because MANE's exon starts at 6534809, but the mapped read starts at 6534806, so the first 3bp of the mapped block fall in MANE's annotated intron and cgranges doesn't return MANE for those 3bp).

#### What does NOT explain the error (disproven)
- Unspliced alignments: ALL mane_absent reads ARE spliced (100% have N in CIGAR) ✗
- Pseudogene contamination: 0 non-GAPDH fragments in GAPDH EM in either oracle or mm2 ✗
- 96.3% tied units: oracle ALSO has 98.9% tied units and quantifies correctly ✗

#### Fix: splice:sr preset + annotation junctions — CONFIRMED EFFECTIVE (2026-03-23)

The original minimap2 command used `-ax splice` (long-read preset, NOT for short reads) and no annotation. The fix is:

```bash
minimap2 -ax splice:sr \
    -j junctions.bed \        # BED12 splice junction annotation
    --secondary=yes -N 20 \
    -t 8 \
    genome.fasta R1.fq.gz R2.fq.gz | \
  samtools sort -n -o output.bam
```

- `splice:sr`: short-read RNA-seq splice preset (added minimap2 v2.29+). Use this for short paired-end RNA-seq, NOT `splice`.
- `-j junctions.bed`: BED12-format annotation file with known splice junctions. Provides anchor boosting at annotated junction sites.
- `junctions.bed` can be generated by `rigel index` (already built for the test case, 40 MB for full human transcriptome).

**Quantitative results comparing all conditions:**

| Transcript | Oracle TPM | OldMM2 (-ax splice) | NewMM2 (splice:sr -j) | Old ratio | New ratio |
|---|---:|---:|---:|---:|---:|
| ENST00000229239.10 (MANE) | 688,435 | 411,208 | **655,241** | 0.597 | **0.952** |
| ENST00000396856.5 | 144,582 | 159,263 | 140,499 | 1.102 | 0.972 |
| ENST00000396861.5 | 77,491 | 58,749 | 69,901 | 0.758 | 0.902 |
| ENST00000396859.5 | 33,921 | 210,821 | **35,951** | 6.215 | **1.060** |
| ENST00000396858.5 | 33,438 | 27,037 | 32,582 | 0.809 | 0.974 |
| ENST00000619601.1 | 22,134 | 108,898 | **29,461** | 4.920 | **1.331** |

The dominant error (+2,659 mane_absent, +1,650 mane_dominated) is almost entirely fixed by using the correct alignment preset and annotation. MANE ratio improves from 0.597 → **0.952** (near perfect). The 6× over-count of ENST396859.5 collapses to 1.06×.

The remaining small residual (0.952 vs 1.0 for MANE; 1.33× vs 1.0 for ENST619601.1) likely reflects reads where the fragment's 5′ end has only 1–3bp in exon 1 before the 240bp intron — too short an anchor even with annotation. This may require the 2-pass alignment approach (write junctions on first pass, boost on second pass) which requires minimap2 v2.31+.

**TODO: Update simulation harness and user documentation to use `splice:sr -j junctions.bed` as the recommended minimap2 command for rigel.**




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

