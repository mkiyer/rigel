# TODO


## Any vestigial / stale hard cutoffs for strand specificity?

Previous versions of the code had hard cutoffs for strand specificity. Do those still exist? These should be removed because our new initialization scheme uses both coverage and strand specificity together. We should not longer require cutoffs for strand specificity anywhere in the code. Do a deep dive thorough search for any remaining strand specificity conditional statements based on some cutoff. 

## gDNA siphoning at imperfect strand specificity

When ss < 1.0, gDNA begins to siphon away fragments. gDNA siphoned ~150/10000 fragments in one case. Why does this happen? With imperfect strand specificity (ss=0.90), mRNA and nRNA may be either sense (90%) or antisense (10%). Both sense and antisense fragments compete for gDNA.

When gDNA contamination is present, the EM **over-assigns** unspliced exonic fragments to gDNA, stealing them from nRNA (and to a lesser extent mRNA). The effect scales with gDNA abundance:

| gDNA level | nRNA rel error | mRNA rel error | gDNA excess |
|-----------|---------------|---------------|-------------|
| 0 | -0.6% | near 0% | 0 |
| 20 | -1.1% | -0.1% | -3.7% |
| 100 | -5.6% | -6.7% | -0.6% |
| 200 | -15.9% | -8.7% | +0.6% |
| 500 | -39.7% | -18.0% | +1.3% |

nRNA is hit harder than mRNA because their fragments overlap more with gDNA (both are unspliced and span intronic regions).

### Root cause: effective length asymmetry

This is **expected behavior**, not a bug. Here's the mechanism:

In the EM, each fragment's posterior probability for component $c$ is:

$$P(c \mid \text{frag}) \propto \frac{\theta_c \cdot P(\text{frag} \mid c)}{\text{eff\_len}_c}$$

For an unspliced exonic fragment (footprint ~200bp) in the diagnostic case (g=500, n=100, m=100):

| Component | profile\_len | eff\_len | $-\log(\text{eff\_len})$ |
|-----------|------------|---------|-------------------------|
| mRNA (exonic) | 4,000 | 3,801 | -8.24 |
| nRNA (genomic span) | 8,000 | 7,801 | -8.96 |
| **gDNA (locus span)** | **17,000** | **16,801** | **-9.73** |

The effective length correction $-\log(\text{eff\_len})$ penalizes components proportionally to their size — a fragment could land anywhere in the component's footprint, so longer components have lower per-position density. **This is correct**: gDNA really does span the whole locus, so a fragment at any position is equally likely.

The siphon happens because:

1. **gDNA has a huge $\theta$**: When g=500, ~89% of fragments are gDNA, so $\theta_{\text{gDNA}} \gg \theta_{\text{nRNA}}$
2. **The eff\_len penalty is sublinear** (logarithmic): gDNA's locus span is only 2× the nRNA span, but its $\theta$ is 10-50× larger. The $\log(\theta)$ advantage overwhelms the $\log(\text{eff\_len})$ disadvantage.
3. **Unspliced exonic fragments are ambiguous**: They're consistent with all three hypotheses (mRNA, nRNA, gDNA). The EM splits them proportional to $\theta / \text{eff\_len}$, which favors gDNA when gDNA dominates.

### Why nRNA is hit harder than mRNA

- **Spliced fragments are safe**: Fragments with annotated splice junctions get $-\infty$ gDNA log-likelihood (gDNA can't produce splice junctions). These are unambiguously mRNA.
- **Intronic fragments are safe**: Fragments overlapping introns are consistent with nRNA and gDNA but NOT mRNA. However, intronic fragments are also ambiguous with gDNA.
- **Exonic unspliced fragments**: Ambiguous across all three. The EM's $\theta$-weighted split siphons some to gDNA. Since nRNA's "claim" on these fragments is weaker (larger eff\_len than mRNA), it loses more.

### Can we minimize it?

This is fundamentally a **model identifiability problem**. When gDNA is very high, the exonic unspliced fragments are genuinely ambiguous — you can't tell from the fragment alone whether it came from nRNA or gDNA. Some approaches to consider:

1. **Already handled well at SS ≥ 0.9**: Strand information resolves most ambiguity. At SS=1.0 and g=100, errors are only ~5-7%. The siphon is really only problematic at very high gDNA (g ≥ 200) relative to RNA.

2. **The nrna\_frac prior helps**: The hierarchical Beta prior anchors the nRNA fraction estimate from density ratios (intronic vs exonic). This already works — the nrna\_frac prior for t0 was 0.633 (close to truth). The issue is that the EM's fragment-level competition overrides this because gDNA fragments vastly outnumber RNA fragments.

3. **Potential improvements** (diminishing returns territory):
   - **Fragment position information**: gDNA fragments should be uniformly distributed across the locus, while nRNA fragments should cluster around gene bodies. Using positional bias models could help disambiguate.
   - **Stronger nrna\_frac prior**: Increasing `kappa_tss` would keep nrna\_frac closer to the density-inferred value, reducing the EM's ability to siphon. But this risks overfitting the prior when it's wrong.
   - **gDNA rate cap**: Bounding gDNA's $\theta$ based on intergenic fragment density could prevent it from absorbing too many genic fragments. This is somewhat already done via the EB gDNA prior.

4. **The practical question**: At g=500 with m=100 and n=100, gDNA is 89% of all fragments. Accurate deconvolution at this contamination level is inherently difficult for any method. The errors at g ≤ 100 (the more realistic range) are modest (< 7%).

**Bottom line**: The gDNA siphon is expected behavior from a well-specified generative model, not a bug. It reflects genuine statistical ambiguity when gDNA dominates. The fix is upstream: reduce gDNA contamination in the wet lab, or use higher strand specificity protocols.



 ## Incorrect nRNA and gDNA initialization

 I agree that there is a bug where the compute_nrna_init() and compute_eb_gdna_priors() are not correctly accessing the strand model. The goal of maintaining multiple strand models (spliced, exon, intron, intergenic) is to facilitate the initialization estimates for seeding the EM. If the spliced strand-specificity (gold standard) is 0.99 and the exonic strand specificity is 0.95 and intronic ss = 0.6, it implies something important about our gDNA contamination rates. We might need to revisit how we estimate gDNA contamination rates (locally, chromosome, globally) at a later time.

## Nascent RNA versus Genomic DNA

In the latest version of the code, nascent RNA appears to outcompete with genomic DNA for fragments. When genomic DNA is high, nascent RNA appears to be siphon off fragments that should be categorized as gDNA.



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


## Fragment length model

- gDNA and RNA should have different models
- different models should be used to help predict
- output shouldn't necessarily report exonic, intronic, intergenic, etc. that's just for initialization



## Output 

- rename 'tpm' -> 'mrna_tpm' and move next to mrna columns
- for tsv files, truncate floating point for readability


## Edit distance

- when reads overlap, edit distance can get double counted in the overlap portion of the reads
- STAR avoids this with the 'nM' tag
- Otherwise, need to add this solution



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