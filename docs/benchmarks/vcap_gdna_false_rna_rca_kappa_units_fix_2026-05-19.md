# VCaP gDNA-to-RNA Error RCA - 2026-05-19

Benchmark output: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/kappa_units_fix`

Truth labels come from flowcell IDs:

- RNA: `C6EL5ANXX`
- gDNA: `H7MFFDSXY`

## Headline

The remaining VCaP error is not primarily a calibration-exposure bug, an aligner
artifact, or a small set of bad coordinates. It is mostly genic, ambiguous,
unspliced gDNA evidence being accepted as RNA evidence.

Post-fix confusion:

| True source | Pred RNA | Pred gDNA | Unresolved |
| --- | ---: | ---: | ---: |
| RNA | 95.48% | 3.64% | 0.89% |
| gDNA | 13.13% | 85.24% | 1.63% |

There are 2,393,558 true-gDNA fragments called RNA. The split is 1,062,195 mRNA
and 1,331,363 nRNA.

## What The Errors Are

Across true-gDNA fragments called RNA:

- 95.95% are `ZS=unspliced`.
- 96.94% are `ZC=ambig_same_strand` or `ZC=ambig_opp_strand`.
- Only ~2.3% are splice-supported (`spliced_annot`, `spliced_implicit`, or
  `spliced_unannot`).
- For mRNA false calls, 96.77% are `all_M`, 98.87% are `NH=1`, and 83.33% have
  `NM=0`.
- For nRNA false calls, 96.56% are `all_M`, 97.79% are `NH=1`, and 81.63% have
  `NM=0`.

So these are mostly clean, uniquely mapped contiguous genomic fragments, not bad
spliced alignments or multimapping noise.

Top error categories:

| Pred detail | ZC | ZS | Count | Fraction |
| --- | --- | --- | ---: | ---: |
| nRNA | ambig_same_strand | unspliced | 877,209 | 36.65% |
| mRNA | ambig_same_strand | unspliced | 830,287 | 34.69% |
| nRNA | ambig_opp_strand | unspliced | 342,908 | 14.33% |
| mRNA | ambig_opp_strand | unspliced | 197,129 | 8.24% |

## Locus Structure

Mega-components matter, but they are not the whole story.

- Global EM locus `3` contributes 408,245 false RNA calls, 17.06% of all
  gDNA-to-RNA errors.
- That locus has 42,639 transcripts, 22,701 nRNA entities, 4,005 genes, and
  2.82M EM fragments.
- Loci with span >= 1 Mb contribute 505,360 false RNA calls, 21.11% of the total.
- Loci with >= 100 transcripts contribute 736,777 false RNA calls, 30.78%.

The error is therefore partly mega-locus RNA-mass borrowing, but most of the
remaining error is distributed across ordinary genic ambiguous loci.

## Local RNA Evidence

False RNA calls are enriched where local RNA evidence exists, but they are not
confined to those windows.

By true RNA-source fragments per 10 kb window:

| RNA-source fragments | gDNA false RNA | gDNA false-RNA rate |
| --- | ---: | ---: |
| 0 | 302,704 | 11.07% |
| 1-10 | 280,725 | 11.51% |
| 11-100 | 417,010 | 12.19% |
| 101-1000 | 883,195 | 15.45% |
| >1000 | 509,924 | 21.69% |

By observed splice-supported RNA calls per 10 kb window:

| Splice-supported RNA calls | gDNA false RNA | gDNA false-RNA rate |
| --- | ---: | ---: |
| 0 | 326,961 | 12.76% |
| 1-10 | 417,553 | 12.55% |
| 11-100 | 615,214 | 14.82% |
| 101-1000 | 759,344 | 14.86% |
| >1000 | 274,486 | 18.16% |

This points to two mechanisms:

1. Local expression pulls some ambiguous genic gDNA fragments into RNA.
2. In windows with little or no splice-supported RNA, ambiguous unspliced
   fragments can still self-seed RNA/nRNA components.

## Guard Experiments

These are post-hoc diagnostics, not production changes. They reassign selected
ambiguous unspliced RNA winners to gDNA and recompute flowcell confusion.

Baseline:

- gDNA -> RNA: 13.13%
- RNA -> gDNA: 3.64%

Naive ambiguous-unspliced guard by winner posterior is too blunt:

| Guard | gDNA -> RNA | RNA -> gDNA |
| --- | ---: | ---: |
| reassign all ambiguous unspliced RNA winners | 0.80% | 29.21% |
| reassign if `ZW < 0.70` | 3.93% | 16.00% |
| reassign if `ZW < 0.50` | 6.37% | 11.67% |

Local splice-support guards are more informative:

| Guard | gDNA -> RNA | RNA -> gDNA |
| --- | ---: | ---: |
| any ambiguous unspliced RNA, no splice-supported calls in 10 kb | 11.46% | 4.24% |
| any ambiguous unspliced RNA, <=10 splice-supported calls in 10 kb | 9.34% | 5.72% |
| nRNA ambiguous unspliced only, no unambiguous splice support in 10 kb | 9.74% | 4.46% |
| nRNA ambiguous unspliced only, <=10 unambiguous splice support in 10 kb | 7.99% | 5.00% |

The nRNA-specific local support guard is the best diagnostic tradeoff in this
round. It directly attacks the larger error class while avoiding the biggest hit
to real mRNA-source fragments.

## Diagnosis

Rigel is still treating ambiguous unspliced genic evidence as too strong for RNA,
especially nRNA. For hybrid-capture data, exonic and exon-adjacent gDNA fragments
are abundant and often cleanly align to transcript-compatible sequence. Without
splice evidence, transcript-unique exonic geometry, or local RNA-only support, a
contiguous genic fragment is biologically gDNA-compatible. The current EM can let
RNA/nRNA mass in the locus convert these fragments into confident RNA winners.

This is most damaging in large connected components, where RNA mass can be
borrowed over a very wide region, but the same failure appears across ordinary
genic loci.

## Recommended Next Target

Target ambiguous unspliced nRNA evidence first.

Concrete next implementation target:

1. Add a local evidence tier for RNA/nRNA calls:
   splice-supported RNA, transcript-unique exonic RNA, and ambiguous-unspliced
   RNA.
2. Downweight or gate ambiguous-unspliced nRNA evidence when a local window or
   partition lacks splice-supported / unambiguous RNA support.
3. Keep mRNA handling more conservative than the naive guards above; mRNA guards
   have a much higher RNA-source cost.
4. Add BAM diagnostics for runner-up gDNA posterior or RNA/gDNA posterior ratio,
   because `ZW` alone cannot distinguish confident RNA-vs-gDNA decisions from
   confident choices within RNA.
5. Separately, regularize or split mega-components so ambiguous unspliced
   fragments cannot borrow RNA abundance from chromosome-scale components.

Priority order:

1. nRNA ambiguous-unspliced local-support prior/gate.
2. Mega-locus localization/regularization.
3. Better diagnostic output for RNA-vs-gDNA posterior margins.

## Artifacts

- Hotspot report: `docs/benchmarks/vcap_gdna_false_rna_hotspots_kappa_units_fix_2026-05-19.md`
- Hotspot tables: `results/vcap_gdna_false_rna_hotspots_kappa_units_fix_2026-05-19/`
- Guard threshold evaluation: `results/vcap_ambiguous_unspliced_guard_eval_2026-05-19/`
- Alignment feature profile: `results/vcap_false_rna_alignment_features_2026-05-19/`
- Local splice support profile: `results/vcap_false_rna_local_splice_support_2026-05-19/`
- Local splice guard evaluation: `results/vcap_local_splice_guard_eval_2026-05-19/`