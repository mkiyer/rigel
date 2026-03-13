# TODO


## Genome partitioning into regions

## Mapping fragments onto regions


The Fragment-Region mapping can be summarized as 'exon_pos_bp', 'exon_neg_bp', 'tx_pos_bp', 'tx_neg_bp'. To compute these sums, we need the Fragment-Region overlap (number of bases) and the Region flags. For example a region (ref='chr1', start=1000, end=2000) with flag (exon_pos = TRUE, exon_neg=FALSE, tx_pos=TRUE, tx_neg=FALSE) overlaps with a Fragment F with aligned block at (ref='chr1', is_reverse=FALSE (pos strand), start=1900, end=2050) overlaps by 100 bases. This can be summed into the Fragment-Region mapping with exon_pos_bp += 100, tx_pos_bp += 100. This should give you an idea of how the Fragment to Region mapping and tabulation/summation should happen. After we compute the Fragment-Region mapping summation/tabulation across all Fragment alignment blocks and their overlapping regions, we then have to determine how to utilize the fragment. When the Fragment-Region mappings are highly ambiguous


## Fragment length model



## Fragment length distribution integration

- this is planned
- strategy to weight observations by uncertainty
- build RNA and gDNA fragment length distributions
- use distributions to improve deconvolution

- gDNA and RNA should have different models
- different models should be used to help predict
- output shouldn't necessarily report exonic, intronic, intergenic, etc. that's just for initialization


Rigel separates gDNA from RNA using splice status and strand specificity. Fragment length is currently unused for discrimination — a single global_model is used for all pools. Bioanalyzer/tapestation confirms gDNA and RNA have distinct insert size distributions.

Key insight: We know strand specificity (ss) and that gDNA is symmetric and perfectly unstranded (0.5). Using Bayes' theorem, we can estimate the gDNA fraction per compartment, then decompose existing fragment length histograms into weighted gDNA and RNA contributions.

Approach: Bayesian Histogram Mixing ("0th-Iteration E-Step")
Instead of arbitrary thresholds or algebraic decontamination, we compute the probability that fragments in each compartment are gDNA, then use those probabilities as fractional weights to build two model histograms via O(1) numpy operations on pre-collected per-compartment histograms.

The Math (per compartment, e.g. "unspliced genic"):
Given S (sense count) and A (antisense count) in a compartment:
R = max(0, (S - A) / (2·s_RNA - 1))     # total RNA
G = (S + A) - R                         # total gDNA
W_sense_gDNA  = (G · 0.5) / S           # gDNA weight for sense hist
W_anti_gDNA   = (G · 0.5) / A           # gDNA weight for antisense hist
Final models:
H_gDNA = H_intergenic + W_sense · H_sense + W_anti · H_anti
H_RNA  = H_spliced    + (1-W_sense) · H_sense + (1-W_anti) · H_anti

Properties:

Zero parameters (no SS thresholds, no min_obs cutoffs)
Graceful degradation: when s_RNA → 0.5, denominator → 0, all W → 0, models collapse to H_intergenic and H_spliced (the only safe anchors)
O(1) mixing: operates on histograms of size 2001, not millions of fragments

Architecture-preserving: frozen before EM, locus-level parallelism untouched

Implementation Steps
Step 1: Collect sense/antisense unspliced histograms during BAM scan

We need 4 raw histograms (2 already exist, 2 new):

H_spliced = category_models[SPLICED_ANNOT] — exists
H_intergenic = intergenic — exists
H_unspliced_same_strand — NEW: unspliced unique-mapper genic fragments where exon_strand == gene_strand
H_unspliced_opp_strand — NEW: unspliced unique-mapper genic where exon_strand != gene_strand

(After strand model finalization, same_strand/opp_strand are mapped to sense/antisense based on p_r1_sense.)
Collection strategy: Add two new vectors to FragLenObservations in the C++ BAM scanner (bam_scanner.cpp). At the point where fragment lengths are recorded for genic fragments (~line 1050), additionally record genomic_footprint into same_strand or opp_strand vectors for unspliced unique-mapper fragments with unambiguous strand.
In Python _replay_fraglen_observations(), route these into two new FragmentLengthModel instances on FragmentLengthModels.
Alternative (no C++ changes): Collect from the buffer after scan via a single vectorized numpy pass using existing buffer columns (splice_type, exon_strand, t_indices, num_hits, genomic_footprint). Fast but adds one O(N) pass.
Files: src/rigel/native/bam_scanner.cpp (add 2 vectors ~15 lines), src/rigel/frag_length_model.py (add 2 models), src/rigel/pipeline.py (replay new observations)
Step 2: Add mix_models() to FragmentLengthModels
New method that performs the Bayesian histogram mixing:
pythondef mix_models(self, s_rna: float, p_r1_sense: float) -> None:
    """Build gDNA and RNA models via Bayesian histogram mixing.

    Uses strand specificity to compute fractional gDNA weights per
    compartment, then mixes sense/antisense unspliced histograms
    with intergenic (pure gDNA) and spliced (pure RNA) anchors.
    """
    # Map same_strand/opp_strand → sense/antisense based on protocol
    if p_r1_sense >= 0.5:
        sense, anti = self.unspliced_same_strand, self.unspliced_opp_strand
    else:
        sense, anti = self.unspliced_opp_strand, self.unspliced_same_strand

    S, A = sense.total_weight, anti.total_weight

    if (2 * s_rna - 1) > 1e-6 and S + A > 0:
        R = max(0.0, (S - A) / (2 * s_rna - 1))
        G = (S + A) - R
        w_s = (G * 0.5 / S) if S > 0 else 0.0
        w_a = (G * 0.5 / A) if A > 0 else 0.0
    else:
        w_s, w_a = 0.0, 0.0  # unstranded fallback

    # gDNA model = intergenic + weighted unspliced
    self.gdna_model.counts = (
        self.intergenic.counts.copy()
        + w_s * sense.counts + w_a * anti.counts
    )
    self.gdna_model._total_weight = float(self.gdna_model.counts.sum())
    self.gdna_model.n_observations = int(self.gdna_model._total_weight)

    # RNA model stored in a new field (or reuse an existing model)
    self.rna_model.counts = (
        self.category_models[SpliceType.SPLICED_ANNOT].counts.copy()
        + (1 - w_s) * sense.counts + (1 - w_a) * anti.counts
    )
    self.rna_model._total_weight = float(self.rna_model.counts.sum())
    self.rna_model.n_observations = int(self.rna_model._total_weight)
Also add rna_model = FragmentLengthModel(max_size=max_size) to __init__.
File: src/rigel/frag_length_model.py
Step 3: Wire mixing into pipeline.py
pythonstrand_models.finalize()
frag_length_models.mix_models(
    s_rna=strand_models.strand_specificity,
    p_r1_sense=strand_models.exonic_spliced.p_r1_sense,
)
frag_length_models.finalize()  # caches LUTs for gdna_model, rna_model, etc.
```

**File**: `src/rigel/pipeline.py` (~line 747-748)

### Step 4: Add gDNA LUT fields to `FragmentScorer`

Add after line 112:
```
gdna_fl_log_prob: np.ndarray | None
gdna_fl_max_size: int
gdna_fl_tail_base: float
File: src/rigel/scoring.py
Step 5: Update from_models() to build both LUTs
Change the RNA LUT source from global_model to rna_model (the mixed RNA model). Add gDNA LUT extraction from gdna_model (the mixed gDNA model). Fall back to global_model if either mixed model has very few observations (< some small threshold like 10).
C++ NativeFragmentScorer automatically picks up the new RNA LUT since fl_log_prob is passed from Python. No C++ constructor changes.
File: src/rigel/scoring.py (from_models(), ~line 162)
Step 6: Add gdna_frag_len_log_lik() function
Same structure as frag_len_log_lik() (line 253) but reads ctx.gdna_fl_* fields.
File: src/rigel/scoring.py (after line 262)
Step 7: Update score_gdna_standalone() and _gdna_log_lik()

score_gdna_standalone() (line 349): use gdna_model.log_likelihood() instead of global_model
scan.py:_gdna_log_lik() (line 429): use gdna_frag_len_log_lik(ctx, ...) instead of frag_len_log_lik(ctx, ...)

Files: src/rigel/scoring.py, src/rigel/scan.py
Step 8: Logging
Enhance frag_length_models.log_summary() and pipeline logging:

Report computed gDNA weights (W_sense, W_anti)
Report gDNA model stats (mean, mode, n_obs) vs RNA model stats
Report whether mixing was applied or fell back (s_rna ≈ 0.5)

Files: src/rigel/frag_length_model.py, src/rigel/pipeline.py
Step 9: Tests
New tests/test_gdna_frag_length.py:

Mixing math: given known S, A, s_rna, verify W_sense and W_anti are correct
Mixing produces distinct gDNA and RNA histograms from synthetic data
Graceful fallback: s_rna ≈ 0.5 → gDNA = intergenic only, RNA = spliced only
gdna_frag_len_log_lik() returns different values from frag_len_log_lik()
score_gdna_standalone() uses gDNA model
End-to-end: pool-specific scoring improves gDNA/RNA discrimination

Update existing scoring tests if they assert specific gDNA/RNA log-likelihood values.
Files Modified
FileChangesrc/rigel/native/bam_scanner.cppAdd 2 new observation vectors (same/opp strand unspliced lengths)src/rigel/frag_length_model.pyAdd rna_model, unspliced_same_strand, unspliced_opp_strand models; add mix_models()src/rigel/scoring.pygDNA LUT fields, gdna_frag_len_log_lik(), update from_models() (RNA LUT from rna_model), update score_gdna_standalone()src/rigel/scan.py_gdna_log_lik() → gdna_frag_len_log_lik()src/rigel/pipeline.pyReplay new observations, call mix_models(), loggingtests/test_gdna_frag_length.pyNew test file
Key Properties

Zero parameters: No SS thresholds, no min_obs cutoffs, no decontamination coefficients
Mathematically closed: The weights balance exactly (expected sense = G·0.5 + R·s_RNA = S)
Graceful degradation: Unstranded → intergenic + spliced anchors only
O(1) mixing: Histogram algebra on arrays of size 2001
Architecture-preserving: Models frozen before EM; locus-level parallelism untouched
All data used: Sense AND antisense unspliced fragments contribute, weighted by their gDNA probability

Verification

pytest tests/ — existing tests pass
pytest tests/test_gdna_frag_length.py — new mixing math + scoring tests
Run on stranded data — log shows distinct gDNA (mean, mode) vs RNA distributions
Compare quantification with/without — no regression, improved gDNA separation



## gDNA siphoning problem

Currently gDNA siphons mRNA + nRNA fragments.

Initialization is a major problem
The EM solver might also need constraints

Two ideas:
1) Predicted expected unspliced vs spliced counts for each transcript. 

We need a function:

unspliced_counts_t <- function(spliced_counts_t, fragment length distribution)

This perfectly complements your Geometric Splicing derivation. The geometric math tells you exactly how many unspliced reads the mRNA *should* have, and these constraints ensure that gDNA can never steal more than its fair, symmetric share of the remainder.

Shall we work out the Python implementation for the Geometric Splicing prediction next?








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



### False positive splice junctions

We can "blacklist" splice junction artifacts by doing genome-wide mapping. Align to the genome and detect cases where gDNA sequence -> spliced alignment (artifact). Compile a "blacklist" from these cases.



## Output 

- rename 'tpm' -> 'mrna_tpm' and move next to mrna columns
- for tsv files, truncate floating point for readability
- nascent rna output




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


## (RESOLVED) Benchmark improvements

- benchmark.py
- salmon needs a salmon index (full transcriptome)
- kallisto needs a kallisto index (full transcriptome)
- build or point to it



## (RESOLVED) nascent RNA decoupling from individual transcripts

In progress

Need to address:
- Output
- Cleanup


Currently each transcript has its own nascent RNA shadow. This is not exactly correct. Transcripts that share the same genomic start and genomic span may have different splicing combinations, but the nascent RNA is by definition the same!

Scenario: 
- Genome 20kb
- Transcript T1 + strand with exons [(1000, 2000), (5000,5500), (7000,7500), (9000,10000)]
- Transcript T2 + strand with exons [(1000, 2000), (9000,10000)]
- Transcript T3 + strand with exons [(1000, 2000), (5000,5500), (9000,10000)]
- Transcript T4 + strand with exons [(4500,5500), (9000,10000)]

Consider this scenario - the current algorithm will define 4 nascent RNA shadows for the 4 transcripts. However, transcripts T1, T2, and T3 actually share the same nascent RNA defined by the transcript span (1000, 10000). The 4 nascent RNAs are exactly the same. The only way nascent RNAs can be distinguished is if their genomic span is different.

This is good news for the Rigel tool! This is because we can consolidate nascent RNAs which will lead to fewer candidates in the EM, and will simplify the EM solver.

We need a complete redesign to support the new definition of nascent RNA. With the new definition, nascent RNA will no longer be associated with each individual transcript. Rather, nascent RNA will be associated with genome_start:genome_end spans.

Please plan the new nascent RNA feature so we can implement it. We need to determine an optimal implementation. Here are my thoughts:

### Implementing the new index

It might be efficient to define nascent RNA at index build time, and have a nascent RNA file (feather, TSV) that can be referenced. The transcript index file will need to be changed as well, to include the nascent RNA index corresponding to each transcript. The index files should be useful to 'rigel' but also downstream tools, so that we can map between transcript <-> nascent RNA <-> gene

### Define nascent RNAs at index time

- At index generation time we can iterate through Transcripts sorted by (ref), genome_start). Nascent RNAs are defined by unique (ref, strand, start, end) tuples in genome coordinate space. 
- We can define the complete list of nascent RNAs at this time and use nrna_id (unique index), ref, strand, start, end, gene index (each gene -> multiple nascent RNAs and each nascent RNA -> 1 gene), perhaps other fields?
- Transcripts will need to be associated with their nascent RNA through a t_to_nrna_index array that maps transcript id -> nascent rna id
- Nascent RNA -> transcripts list? Not sure if we need to maintain this or not but it would be key by nrna_idx -> value is list of transcript_ids.. i have to think whether the algorithm will need this

### Change nascent RNA initialization

Currently to initialize nascent RNA shadows there must be evidence of fragment overlapping an "unambiguous intron" interval. In otherwords, a fragment must at least partially overlap an intronic interval where there is no exon on either strand. This is sufficient to initialize nascent RNA shadows.

### "linking" mature RNA to nascent RNA

Currently we model the steady state nascent fraction for each transcript using a nrna fraction 'beta'. This nascent RNA fraction  remains a 'transcript level' property. 

The derived abundance for mature RNA is the same: 
mature RNA = theta_t (from EM) x (1 - nrna_fraction_t)

However, the abundance for each nascent RNA different:

nascent RNA = sum((theta_t x nrna_fraction_t) for t in (transcripts that map to this nascent RNA))

This models nascent RNA as a "one-to-many" definition: one nascent RNA can produce multiple mature RNAs depending on the splicing fraction. 

### Implementation

This is a major architectural change. It needs to be implemented carefully, in stages. We must assure no regressions at each stage. Ultimately this is a more biologically correct definition of nascent RNA and should lead to more efficient EM behavior (we make the EM assign nascent RNA to transcripts)

### Backward compatibility

There is NO need to support the prior code base or method. We can forge ahead, delete/discard older methods, and produce production quality CLEAN code without stale/vestigial methods or redundancies.


Now, consider this feature and PLAN the implementation in phases. Write an extremely detailed implementation plan in 'docs' named 'nascent_rna_decoupling.md' 


### EM algorithm

The number of components in the EM algorithm will change:

[0,T) - mature RNA components (stays the same)
[T,T+N) - nascent RNA components (where 'N' is number of nascent RNAs)
T+N - genomic DNA (stays the same, single component)

EM M-Step -- The splicing fraction parameters are updated during the 'M' step of the EM algorithm. The splicing fraction parameters remain specific to each transcript, that doesn't change.


### Output

Currently we output transcript abundances (separate file), gene abundances estimated from the transcripts (separate file). Nascent RNA levels are included in the transcript file. Now, nascent RNA is no longer defined as an individual transcript property. It is associated with a group of transcripts that share the same genomic start/end coordinates.

Therefore, we will need to create a new nascent RNA output file. It can include: nrna_idx (index), gene_id, gene_name, transcript_ids (not sure whether to include, it would be a comma-delimited list, but the transcript file could be used to lookup the nascent RNA without this), effective length, counts, and abundance (tpm)

The transcript abundance output file will change. It will no longer report nascent RNA counts. It should add: nrna_idx (index of nascent RNA that produces this mature RNA), nascent rna fraction (0.0-1.0) fraction of nascent RNA.

The gene-level output file can be updated as well. In the gene file, we should add gene-level nascent estimation which can be a weighted average of the nascent RNA counts (same computation as we do for transcripts). This should be added.
