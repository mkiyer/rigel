# TODO


## BAM scan pileup

I have an idea that is an offshoot of your "5. Coverage-shape coherence priors" -- during the 

## gDNA -> RNA failures



## increasing region granularity

One of Rigel's main goals is to correctly deconvolute genomic DNA (gDNA, unspliced, unstranded), nascent RNA (nRNA, unspliced, strand-specific), and mature RNA (mRNA, spliced and unspliced, strand-specific).

We continue to have higher than acceptable error. This PR attempts to decrease pool-level error, correctly resolving gDNA, nRNA, and mRNA.

The tool partitions the genome into non-overlapping regions. The current behavior is to cluster overlapping exons into a single conglomerated region. This mostly works well, but occasionally yields some large exonic regions that contain multiple overlapping transcripts and genes. Sometimes the culprit is a very long single-exon transcript. Other times it is can be read-through transcripts that splice from one gene into another. Whatever the cause, large exons destroy granularity that we need. We hypothesize that restoring granular regions (instead of coarse regions) will improve the tool behavior overall.

This PR is to refactor our region partioning and the downstream analysis. This PR will touch many parts of the code including the Rigel index, BAM scanner/buffer, scoring, calibration, locus building, and setting the bayesian prior and warm start for EM.

Fortunately, Rigel USED to have a granular region system in place but we abandoned this in favor of the current coarse but simpler region system.

You can audit and retrieve code for the prior region system from the git repo as commit '14c307f'.

Our new granular region system, which I want you to reimplement, is defined as follows:

A region is defined as a genomic interval with (ref_id, start, end). We will partition the entire genome into non-overlapping regions.

Each region is defined by a set of true/false flags:
- intron_pos: an intron on the genomic positive strand overlaps this region
- intron_neg: an intron on the genomic negative strand overlaps this region
- exon_pos: an exon on the genomic pos strand overlaps this region
- exon_neg: an exon on the genomic neg strand overlaps this region

This 4-bit, 4-flag signature defines a region. Note: we may have used tx_pos and tx_neg instead of intron_pos, intron_neg in the past. These are interchangable; one can be computed from the other.

if bits are intron_pos, intron_neg, exon_pos, exon_neg

b0000 = intergenic
b0001 = exon neg
b0010 = exon pos
b0011 = overlapping exons on pos and neg strands
b0100 = intron pos
b0101 = intron on pos strand overlapping with an exon on the neg strand
b0110 = intron on pos strand overlapping exon on pos strand

.... and so on and so forth hopefully you understand the pattern here?
All combinations of 4 bits are possible


The region partitioning can then be defined. Given a genome (fasta file) and gene annotation (gtf file):

- For each reference chromosome:
   - create a sorted list of unique exon boundaries by iterating over transcript exons. each boundary has the following info: (is_tss = is the boundary of a transcript start site, is_tes = is the boundary of a transcript end site, is_exon_start, is_exon_end)
   - initialize 4-bit state as b0000 (intergenic) and iterate from left to right over sorted boundaries
      - if state *changes*, then save current region and its state; update 4-bit state for each boundary encountered

This is roughly how the region partition definition algorithm works.

Today's "coarse" region partitioning can be easily derived from the fine-grained region partitioning.

The coarse partioning is unstranded and only has INTERGENIC, INTRON, and EXON classes. We also have an ambiguous (AMBIG) class.

if bits are 3:intron_pos, 2:intron_neg, 1:exon_pos, 0:exon_neg

4-bit states new fine-grained system -> old coarse system
b0000 = INTERGENIC
b0001 = EXON
b0010 = EXON
b0011 = EXON (AMBIG)
b0100 = INTRON
b0101 = EXON
b0110 = EXON (AMBIG)
b0111 = EXON (AMBIG)
b1000 = INTRON
b1001 = EXON (AMBIG)
b1010 = EXON
b1011 = EXON (AMBIG)
b1100 = INTRON (AMBIG)
b1101 = EXON (AMBIG)
b1110 = EXON (AMBIG)
b1111 = EXON (AMBIG)

You'll need to verify and derive this for yourself. This is meant to give you a pattern. Just looking at this table should make it clear why the coarse-grained region partioning obfuscates a lot of important granularity. 

We need to revert our code base to the fine-grained region system. This will give us finer regions to evaluate gDNA density. It should yield better estimates of gDNA density. It should yield better estimates of exposure weight for gDNA effective length normalization.

I need you to do the following:

1) find code in the older git repo that contained the prior region calibration system
2) plan a migration from coarse-grained regions to fine-grained regions. This will need to happen in phases.

Here are proposed phases:

Phase 0: code cleanup. Audit the rigel index, region partioning, bam scanning, and calibration code. Try to clean up aspects that are stale, overly complex, legacy, dead, etc. We want to implement over a clean slate.

**NOTE** We are going to implement this with ZERO backwards compatibility and ZERO legacy compatibility

Phase 1: update rigel index to use fine grained regions. rewrite the rigel index code to derive fine-grained regions instead of coarse-grained regions.

Phase 2: BAM scanning. In the BAM scanner, we accumulate fragments over matching regions. We are going to revamp the way we do this accumulation. The new accumulation will give us powerful new information that will be used during calibration.

When we are scanning the BAM file, we construct Fragment objects. A fragment object is a collection of aligned blocks (ref_id, start, end) and other information. We need to map fragments onto regions. A lot of the code to do this can be reused.

Currently, we exclude multi-mapping, chimeric, splice artifacts, and perhaps other classes of transcripts that are difficult to interpret during the calibration step. Likely we can use the existing fragment exclusion criteria.

For region accumulation, each regions needs to store fragment count data as fractional counts (float). Fragments are either CONTAINED (completely contained within a single region) or CROSSING (overlapping multiple regions, e.g. generates boundary flux)

- CONTAINED unspliced pos strand, unspliced neg strand, spliced pos strand, spliced neg strand
- BOUNDARY FLUX LEFT unspliced pos, unspliced neg, spliced pos, spliced neg
- BOUNDARY FLUX RIGHT unspliced pos, unspliced neg, spliced pos, spliced neg


For eligible fragments, we map them over the region partioning using the compatible aligned (ref_id, start, end) blocks.

- For each fragment (spliced and unspliced):
   - for each aligned block (ref_id, start, end):
      - determine the set of overlapping regions
      - if overlaps a single region update region data for contained
      - else for each overlapping region R_i (region_ref_id, region_start, region_end):
         - compute number of overlapping bases (overlap_bp)
         - fractional overlap is overlap_bp / aligned_block_size (end - start)
         - determine whether left, right, or both boundaries are crossed
         - update boundary flux for each boundary crossed. fractional overlap must be divided by number of boundaries crossed (no double counting)

We required that each fragment produces exactly 1.0 count units that are partioned among contained and boundary crossing compartments. For example, take regions R1 =(0, 100), R2=(100, 150), R3=(150, 300).

Fragment has aligned block (80, 140). Aligned block is 60bp Fragment overlaps R1 and updates R1 boundary right with (20/60)=0.33, then R2 boundary left with (40/60) = 0.66. 

Another example fragment with aligned block (50,250) overlaps R1 right boundary (50/200) = 0.25, R2 both left and right boundaries so R2 left (50/200/2) = 0.125 and R2 right = (50/200/2) = 0.125, then R3 left boundary (100/200) = 0.5. The sum of fragment contributions = 1.0.

This accumulation paradigm accomplishes the following:
- fine-grained
- maintains strand information
- maintains spliced/unspliced information
- allows computation of boundary flux (compatible with today's behavior)
- allows regions to compute their "independent" fragment support by aggregation boundary crossing and contained fragments

I believe that this region accumulation system will support a more sophisticated calibration step.

We can recover existing behavior from the region data if we choose to. But we can also extend existing behavior. We need to develop enhanced calibration that uses fine-grained regions and boundary flux.

Phase 3: Fragment length model

During BAM scanning, we must continue to identify gDNA-compatible fragments for the gDNA FL model. We can keep this approach the same as the current approach (INTERGENIC, INTRONIC, and EXON-INTRON fragments are used for the gDNA model).

We may choose to revert to an older gDNA FL model when we get here.

Phase 4: Calibration

This is not designed yet.


Begin designing this new implementation. Store the design in 'docs/fineregions'. 



## argument for fractional accumulator

I want to discuss fractional 12-channel accumulator that you felt that we should defer indefinitely. Looking through your rationale now -- 1) It is wonderful that the current CalibrationAccumulator is sophisticated. We would KEEP most of the current infrastructure when we change to float32 (not float64) fractional counts. We would expand boundary flux with orientation to include spliced fragments.  2) We would KEEP splicing anchor tolerance (excluding fragments from flux when overlap is less than `tolerance` in bp). The current splicing_anchor_tolerance support is wonderful and should work. 3) Integer accumulator creates a "double counting" problem -- suppose a fragment overlaps 7 fine-grained regions. It would increment boundary flux in all regions, but there would be no way to know that the incremented boundary flux came from a single fragment. We have fragments of many different lengths: an 800bp fragment might overlap many more regions than a 150bp fragment. This creates bias in integer accumulation to longer fragments. For coarse-grained accumulation, we did not pay attention to this but with fine-grained accumulation, we are breaking apart larger EXON regions into many smaller regions. Fragments will often overlap many regions. Why do you feel so strongly about deferring our proposed fractional 12-channel accumulator? This concept would allow each region to independently accumulate fragment support. One of the nice properties is that for a single region, we can sum the boundary flux and contained fractional support in a region and easily obtain a total. The sum across regions will equal the total fragments used in accumulation. Without this, we cannot easily solve "how many fragments are in this region?". Right? ---> In summary, I don't see or understand how an expanded fractional accumulation creates major problems. This would add support for spliced fragment boundary flux (new feature) AND make boundary flux universally helpful rather than just used for EXON-INTRON boundaries. We would use boundary flux everywhere, for all regions. -----> I have not developed the algorithms yet, but I believe that this expanded fractional accumulator formulation will pave the way for more sophisticated gDNA estimation. I am thinking about complex sets of regions with overlapping introns and exons and both pos/neg strands. Instead of one coarse EXON region, we will have many adjacent exonic regions. My plan is to devise a gDNA density estimator for these internal exon regions that incorporates spliced/unspliced boundary flux and contained fragments to predict gDNA densities. I don't think this would be possible without this formulation. If you have a better formulation or can suggest improvements please do! Also, I would love you to do make the implementation extremely efficient! We could certainly boost performance if we think about implementation ideas. This is a first pass at giving us the foundation we need.

### splicing anchor tolerance goes away

Splicing anchor tolerance logic needs to be examined in the setting of our new fractional accumulator, because it may no longer be needed. The splicing anchor tolerance was intended for fragments that cross exon-intron boundaries (exon region -> intron region) but the intronic overlap was very small (<3 bp default). this small "overhang" could be aligner error where that last few bases may truly belong on the next exon (spliced) but the aligner didn't align it correctly. when we are computing exon-intron boundary flux, we want to use that as a proxy for gDNA and the fragments with tiny overhang into the intron would otherwise contaminate our exon-intron signal with RNA. However, let's rethink this because the new fractional accumulation paradigm may elegantly solve the splicing overhang problem without the need for a tolerance. In the new fractional accumulator, a fragment that overlaps an exon-intron region pair and overhangs into an intron would contribute to boundary flux, but that flux would get partitioned between the exon and intron. If the aligned block is 100bp with 98bp exon and 2bp intron overhang, the flux would be partitioned so that 98/100 would be allocated to the exon and 2/100 would be allocated to the intron. So if we use the intron region boundary flux to compute gdna density, we would only see 2/100 contribution from the overhanging read -- it is true that this is rna noise, but far less than the prior case where every boundary crossing fragment gets equal integer weight of 1. The intron-to-exon flux which will be useful for gdna estimation is weighted by intronic overlap. So it may be the case the the splicing anchor tolerance isssue is naturally handled by this approach. Can you evaluate this carefully? Derive the logic yourself and work through some examples. Can we rely upon our fractional accumulator to handle splicing anchor tolerance naturally? If so, this simplifies our algorithm. Go ahead and rework the design document v3 to remove splicing anchor tolerance logic, unless you believe that it still serves an important purpose

## hybrid capture simulator

We have an extensive code base for simulation in place in the 'src/sim' folder within Rigel. To support hybrid capture simulation, the simulator needs additional inputs. 1) It needs to accept a BED file with capture targets. It needs to support simple single block BED format column 1 = ref, column 2 = start, column 3 = end, OR the BED12 format in case we have exon-exon spanning probes (splice junction spanning probes) in which case a single probe has multiple genomic capture intervals. 2) We need to specify an overall enrichment ratio (e.g. 1000X) that conveys how much enrichment the capture probes provide. This is the global on-target versus off-target probability ratio. 3) We need to define the probability landscape for sampling reads. Instead of sampling reads uniformly, we need to sample them according to a genome


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

