# VCaP AR-Locus Subset Analysis - 2026-05-20

Region: `chrX:66929868-69019700`

Source BAM: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam`

Output directory: `/Users/mkiyer/Downloads/rigel_runs/vcap_ar_locus/`

## BAM Outputs

| File | Purpose |
| --- | --- |
| `ar_locus.input.namesort.bam` | Name-sorted AR-locus subset used as Rigel input. |
| `ar_locus.source_annotated.coord.bam` | Coordinate-sorted/indexed subset carrying the source BAM's prior annotation tags. |
| `ar_locus.rigel.annotated.coord.bam` | Coordinate-sorted/indexed annotated BAM from the default sampled assignment run. |
| `ar_locus.rigel.map.annotated.coord.bam` | Coordinate-sorted/indexed annotated BAM from deterministic MAP assignment. Best file for browser review. |
| `ar_locus.rigel.map.gdna_false_rna.coord.bam` | Coordinate-sorted/indexed BAM containing only true-gDNA fragments assigned to RNA in MAP mode. |

## Local Rigel Run

- Input subset: `133,290` alignment records, `64,780` representative primary fragments.
- Truth by flowcell: `44,005` RNA-source fragments, `20,775` gDNA-source fragments.
- Rigel local run completed successfully and emitted quant tables, locus stats, summaries, and annotated BAMs.
- Local AR locus is locus `51`: `13` transcripts, `1` gene, `59,659` EM fragments.

## Confusion

Default sampled assignment:

| True source | Pred RNA | Pred gDNA | Pred unresolved |
| --- | ---: | ---: | ---: |
| gDNA | 4,014 / 20,775 (19.32%) | 16,677 / 20,775 (80.27%) | 84 / 20,775 (0.40%) |
| RNA | 42,683 / 44,005 (96.996%) | 1,267 / 44,005 (2.88%) | 55 / 44,005 (0.12%) |

Deterministic MAP assignment:

| True source | Pred RNA | Pred gDNA | Pred unresolved |
| --- | ---: | ---: | ---: |
| gDNA | 3,931 / 20,775 (18.92%) | 16,760 / 20,775 (80.67%) | 84 / 20,775 (0.40%) |
| RNA | 42,831 / 44,005 (97.33%) | 1,119 / 44,005 (2.54%) | 55 / 44,005 (0.12%) |

The MAP result is almost as bad as the sampled assignment result, so this is a real model preference/ambiguity issue, not mainly sampling noise.

## False-RNA Anatomy In MAP Mode

- `3,649 / 3,931` false-RNA fragments are called mRNA; only `282 / 3,931` are called nRNA.
- `3,372 / 3,931` false-RNA fragments go to `mRNA AR ENST00000396043.4`.
- `3,587 / 3,931` false-RNA fragments are pure `EXON`, dominant `EXON`, `POS` strand.
- The pure EXON/POS class has a `47.50%` gDNA-to-RNA false rate: `3,587 / 7,551`.
- The dominant false class is `mrna / ambig_same_strand / unspliced / NH=1`: `3,612` fragments.
- Top 10 kb false window: `chrX:67,539,868-67,549,867`, with `2,788 / 6,585` true-gDNA fragments assigned to RNA (`42.34%`).

## Root Cause

The AR-locus failure is primarily mature-exon gDNA being confused with highly expressed AR mRNA. These are not mostly splice artifacts, multimappers, or nRNA siphons. They are clean, unique, unspliced, same-strand fragments landing inside AR exons, where mRNA and gDNA have very similar fragment-level evidence.

The local subset also exposes a calibration limitation: Rigel still computes gDNA density against whole-genome opportunity even though the BAM was region-extracted. That makes the local prior extremely weak. In the local AR run, locus `51` has `gdna_prior_count_em=1.38` for `59,659` EM fragments. The EM does recover some gDNA (`14,196` in the sampled run, `16,673` in MAP), but the prior is far too small to prevent exon-contained same-strand fragments from being absorbed by the very abundant AR mRNA component.

This means the next target should not be another nRNA-specific fix. The immediate target is capture/local-context-aware gDNA modeling for exonic regions: the gDNA opportunity denominator and/or regional exposure field needs to know the assay's accessible/captured footprint, so exon-contained gDNA fragments in enriched regions receive a realistic local gDNA prior.