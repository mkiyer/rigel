# VCaP gDNA-to-RNA Error By Region - 2026-05-19

BAM: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/kappa_units_fix/annotated.bam`

Index: `/Users/mkiyer/Downloads/rigel_runs/refs/rigel_index`

This analysis classifies each counted fragment by overlap with Rigel's calibration
`regions.feather` partition. It uses both mates' primary aligned blocks, applies
minimum region overlap `q=3` bp, and then counts the hard Rigel
winner from the fragment-level `ZF` tag.

Fragments counted: 32,218,865

True gDNA fragments: 18,228,677

True gDNA -> RNA calls: 2,393,558 (13.13%)

## Region Headline

- EXON-containing masks explain 2,333,137 gDNA->RNA calls (97.48% of false RNA).
- INTRON-only explains 58,231 calls (2.43%).
- INTERGENIC-only explains 1,829 calls (0.08%).

## Exact Region Mask

| Region mask | Broad class | gDNA total | gDNA->RNA | False rate | False RNA frac |
| --- | --- | --- | --- | --- | --- |
| EXON_ONLY | EXON_CONTAINING | 5,385,581 | 1,459,441 | 27.10% | 60.97% |
| EXON_INTRON | EXON_CONTAINING | 11,130,820 | 871,589 | 7.83% | 36.41% |
| INTRON_ONLY | INTRON_ONLY | 885,205 | 58,231 | 6.58% | 2.43% |
| INTERGENIC_ONLY | INTERGENIC_ONLY | 295,893 | 1,829 | 0.62% | 0.08% |
| EXON_INTERGENIC | EXON_CONTAINING | 264,118 | 1,344 | 0.51% | 0.06% |
| EXON_INTRON_INTERGENIC | EXON_CONTAINING | 115,724 | 763 | 0.66% | 0.03% |
| INTRON_INTERGENIC | INTRON_ONLY | 3,478 | 361 | 10.38% | 0.02% |
| NONE | NONE | 147,858 | 0 | 0.00% | 0.00% |

## Top Dominant Regions

Dominant region means the qualified region with the largest aligned-block overlap
for the fragment. Boundary is `left/right` from `regions.feather`.

| Region | Type | Len bp | Strand | gDNA->RNA | Rate | mRNA | nRNA | Exon tx | Exon genes | Span genes | Boundary |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr2:178,559,311-178,576,831 | EXON | 17,521 | AMBIG | 4,607 | 67.37% | 341 | 4,266 | 16 | 2 | 2 | 1/1 |
| chr11:65,497,640-65,508,073 | EXON | 10,434 | AMBIG | 4,396 | 80.76% | 2,157 | 2,239 | 70 | 5 | 5 | 0/0 |
| chr1:152,348,735-152,357,647 | EXON | 8,913 | NEG | 3,067 | 67.58% | 1,449 | 1,618 | 1 | 1 | 2 | 1/1 |
| chrX:67,544,021-67,546,762 | EXON | 2,742 | POS | 2,793 | 45.25% | 2,711 | 82 | 8 | 1 | 1 | 0/1 |
| chr9:106,924,133-106,929,759 | EXON | 5,627 | POS | 2,204 | 49.68% | 1,419 | 785 | 6 | 1 | 1 | 1/1 |
| chr9:76,704,761-76,711,358 | EXON | 6,598 | AMBIG | 1,860 | 79.08% | 1,033 | 827 | 11 | 2 | 2 | 1/1 |
| chr21:33,549,476-33,555,391 | EXON | 5,916 | POS | 1,860 | 47.12% | 1,851 | 9 | 7 | 2 | 2 | 1/1 |
| chr1:152,212,076-152,221,490 | EXON | 9,415 | NEG | 1,789 | 65.99% | 823 | 966 | 1 | 1 | 2 | 1/1 |
| chr9:32,629,454-32,635,669 | EXON | 6,216 | AMBIG | 1,771 | 49.90% | 240 | 1,531 | 2 | 2 | 2 | 0/1 |
| chr1:152,302,165-152,314,957 | EXON | 12,793 | AMBIG | 1,672 | 61.47% | 762 | 910 | 5 | 2 | 2 | 1/1 |
| chr2:185,788,643-185,797,526 | EXON | 8,884 | AMBIG | 1,619 | 66.03% | 802 | 817 | 6 | 2 | 2 | 1/1 |
| chr11:64,116,219-64,119,210 | EXON | 2,992 | AMBIG | 1,585 | 39.75% | 1,013 | 572 | 3 | 2 | 3 | 1/1 |
| chr2:178,530,241-178,535,849 | EXON | 5,609 | AMBIG | 1,529 | 65.74% | 114 | 1,415 | 40 | 2 | 2 | 1/1 |
| chr11:62,515,902-62,534,074 | EXON | 18,173 | NEG | 1,433 | 24.98% | 1,433 | 0 | 4 | 1 | 1 | 1/1 |
| chr13:23,328,826-23,341,690 | EXON | 12,865 | NEG | 1,398 | 35.06% | 1,398 | 0 | 12 | 1 | 1 | 1/1 |
| chr2:178,744,403-178,752,039 | EXON | 7,637 | AMBIG | 1,312 | 60.49% | 0 | 1,312 | 3 | 2 | 2 | 1/1 |
| chr13:109,782,042-109,786,583 | EXON | 4,542 | AMBIG | 1,268 | 33.12% | 1,260 | 8 | 2 | 2 | 2 | 1/1 |
| chr8:76,850,886-76,856,300 | EXON | 5,415 | POS | 1,182 | 41.86% | 1,135 | 47 | 6 | 1 | 1 | 1/1 |
| chr15:99,128,832-99,135,593 | EXON | 6,762 | AMBIG | 1,176 | 57.06% | 735 | 441 | 5 | 2 | 2 | 1/0 |
| chr8:105,801,047-105,804,539 | EXON | 3,493 | POS | 1,119 | 67.45% | 631 | 488 | 4 | 1 | 2 | 1/1 |
| chr13:102,729,367-102,750,388 | EXON | 21,022 | NEG | 1,107 | 10.47% | 1,107 | 0 | 1 | 1 | 1 | 0/1 |
| chr19:56,810,077-56,817,838 | EXON | 7,762 | AMBIG | 1,105 | 50.00% | 608 | 497 | 23 | 3 | 3 | 1/1 |
| chr7:100,951,841-100,960,645 | EXON | 8,805 | POS | 1,097 | 34.86% | 530 | 567 | 3 | 1 | 1 | 1/1 |
| chr8:143,915,153-143,922,395 | EXON | 7,243 | NEG | 1,091 | 39.88% | 1,074 | 17 | 14 | 1 | 1 | 0/1 |
| chr2:185,799,697-185,809,133 | EXON | 9,437 | AMBIG | 1,058 | 42.35% | 978 | 80 | 5 | 2 | 2 | 1/1 |
| chr5:180,847,611-180,851,493 | EXON | 3,883 | NEG | 1,042 | 37.73% | 1,016 | 26 | 7 | 1 | 1 | 0/1 |
| chr19:8,945,497-8,967,189 | EXON | 21,693 | NEG | 1,031 | 16.38% | 1,031 | 0 | 4 | 1 | 1 | 1/1 |
| chr3:49,716,829-49,722,147 | EXON | 5,319 | AMBIG | 1,010 | 60.84% | 330 | 680 | 17 | 3 | 3 | 1/1 |
| chr10:128,102,579-128,109,423 | EXON | 6,845 | NEG | 990 | 42.38% | 990 | 0 | 3 | 1 | 1 | 1/1 |
| chr11:63,718,702-63,721,032 | EXON | 2,331 | POS | 987 | 29.52% | 987 | 0 | 5 | 1 | 1 | 1/1 |
| chr8:143,857,324-143,873,298 | EXON | 15,975 | NEG | 958 | 25.55% | 116 | 842 | 2 | 1 | 1 | 0/1 |
| chr4:87,613,785-87,616,873 | EXON | 3,089 | POS | 948 | 54.64% | 549 | 399 | 1 | 1 | 2 | 1/1 |
| chr21:14,961,235-14,968,526 | EXON | 7,292 | AMBIG | 944 | 40.03% | 944 | 0 | 6 | 2 | 3 | 1/1 |
| chr1:15,928,091-15,936,266 | EXON | 8,176 | POS | 940 | 39.18% | 903 | 37 | 4 | 1 | 1 | 1/1 |
| chr4:143,695,491-143,700,675 | EXON | 5,185 | AMBIG | 922 | 47.75% | 612 | 310 | 2 | 2 | 2 | 1/1 |
| chr11:66,857,064-66,859,093 | EXON | 2,030 | AMBIG | 920 | 62.59% | 550 | 370 | 4 | 2 | 2 | 1/1 |
| chr11:108,505,435-108,514,875 | EXON | 9,441 | NEG | 873 | 33.11% | 873 | 0 | 5 | 1 | 1 | 0/1 |
| chr2:178,547,407-178,550,681 | EXON | 3,275 | AMBIG | 873 | 69.23% | 37 | 836 | 13 | 2 | 3 | 1/1 |
| chr22:42,209,651-42,215,440 | EXON | 5,790 | NEG | 857 | 39.95% | 857 | 0 | 7 | 1 | 1 | 1/1 |
| chr4:76,738,761-76,742,130 | EXON | 3,370 | POS | 849 | 52.02% | 646 | 203 | 5 | 1 | 2 | 1/1 |

Top dominant-region concentration:

- Top 20 regions: 39,640 (1.66% of false RNA)
- Top 100 regions: 104,203 (4.35% of false RNA)

## Dominant Region Bins

| Type | Len bin | Strand | Boundary | Regions | gDNA total | gDNA->RNA | False rate |
| --- | --- | --- | --- | --- | --- | --- | --- |
| EXON | 1k-5k | POS | 1/1 | 6,794 | 1,070,675 | 245,775 | 22.96% |
| EXON | 1k-5k | NEG | 1/1 | 6,298 | 966,074 | 220,216 | 22.79% |
| EXON | 1k-5k | AMBIG | 1/1 | 3,129 | 581,791 | 165,155 | 28.39% |
| EXON | 101-250 | POS | 1/1 | 42,598 | 1,117,218 | 124,656 | 11.16% |
| EXON | 501-1k | POS | 1/1 | 6,749 | 617,958 | 116,397 | 18.84% |
| EXON | 101-250 | NEG | 1/1 | 42,174 | 1,092,370 | 115,608 | 10.58% |
| EXON | 1k-5k | POS | 1/0 | 3,526 | 519,384 | 105,119 | 20.24% |
| EXON | 501-1k | NEG | 1/1 | 6,657 | 576,842 | 104,753 | 18.16% |
| EXON | 1k-5k | NEG | 0/1 | 3,459 | 500,347 | 102,192 | 20.42% |
| EXON | 251-500 | NEG | 1/1 | 11,461 | 664,166 | 94,171 | 14.18% |
| EXON | 251-500 | POS | 1/1 | 11,274 | 649,962 | 88,490 | 13.61% |
| INTRON | 1k-5k | POS | 0/0 | 41,264 | 1,348,066 | 70,526 | 5.23% |
| INTRON | 1k-5k | NEG | 0/0 | 39,746 | 1,199,217 | 56,741 | 4.73% |
| EXON | 5k-10k | AMBIG | 1/1 | 491 | 162,661 | 51,752 | 31.82% |
| EXON | 5k-10k | POS | 1/1 | 387 | 127,095 | 35,651 | 28.05% |
| INTRON | 101-250 | POS | 0/0 | 9,328 | 292,268 | 32,484 | 11.11% |
| EXON | 5k-10k | NEG | 1/1 | 348 | 110,664 | 29,805 | 26.93% |
| INTRON | 101-250 | NEG | 0/0 | 9,124 | 273,221 | 28,946 | 10.59% |
| EXON | 5k-10k | NEG | 0/1 | 518 | 137,449 | 28,897 | 21.02% |
| INTRON | 501-1k | POS | 0/0 | 13,848 | 483,423 | 28,705 | 5.94% |
| EXON | 501-1k | AMBIG | 1/1 | 1,419 | 117,179 | 27,506 | 23.47% |
| EXON | 5k-10k | POS | 1/0 | 507 | 120,313 | 27,113 | 22.54% |
| INTRON | 501-1k | NEG | 0/0 | 13,982 | 442,821 | 24,711 | 5.58% |
| INTRON | 1k-5k | AMBIG | 0/0 | 9,823 | 241,671 | 24,543 | 10.16% |
| INTRON | 251-500 | POS | 0/0 | 9,845 | 348,350 | 24,250 | 6.96% |
| INTRON | 251-500 | NEG | 0/0 | 9,425 | 308,401 | 19,950 | 6.47% |
| INTRON | 5k-10k | POS | 0/0 | 12,631 | 346,556 | 17,890 | 5.16% |
| INTRON | 10k-50k | POS | 0/0 | 15,113 | 351,299 | 17,762 | 5.06% |
| EXON | 501-1k | POS | 1/0 | 1,215 | 95,222 | 15,889 | 16.69% |
| EXON | 10k-50k | AMBIG | 1/1 | 86 | 47,969 | 15,753 | 32.84% |
| INTRON | 10k-50k | NEG | 0/0 | 14,163 | 322,840 | 15,619 | 4.84% |
| EXON | 501-1k | NEG | 0/1 | 1,162 | 87,531 | 14,561 | 16.64% |
| EXON | 1k-5k | AMBIG | 1/0 | 368 | 54,670 | 14,339 | 26.23% |
| INTRON | 101-250 | AMBIG | 0/0 | 2,525 | 61,945 | 13,983 | 22.57% |
| EXON | 1k-5k | AMBIG | 0/1 | 377 | 57,458 | 13,837 | 24.08% |
| INTRON | 5k-10k | NEG | 0/0 | 11,973 | 302,708 | 13,825 | 4.57% |
| EXON | 1-100 | POS | 1/1 | 16,293 | 110,224 | 13,820 | 12.54% |
| EXON | 1k-5k | POS | 0/1 | 423 | 68,359 | 13,441 | 19.66% |
| EXON | 1k-5k | NEG | 1/0 | 410 | 70,503 | 12,454 | 17.66% |
| EXON | 1-100 | NEG | 1/1 | 15,419 | 99,831 | 12,309 | 12.33% |
| INTRON | 501-1k | AMBIG | 0/0 | 3,826 | 100,848 | 11,840 | 11.74% |
| INTRON | 251-500 | AMBIG | 0/0 | 2,801 | 78,888 | 11,145 | 14.13% |
| EXON | 251-500 | AMBIG | 1/1 | 1,116 | 60,964 | 10,608 | 17.40% |
| EXON | 1k-5k | NEG | 0/0 | 752 | 89,455 | 10,098 | 11.29% |
| EXON | 1k-5k | POS | 0/0 | 774 | 83,233 | 8,453 | 10.16% |
| EXON | 501-1k | NEG | 1/0 | 819 | 61,297 | 8,233 | 13.43% |
| EXON | 501-1k | POS | 0/1 | 847 | 57,288 | 7,812 | 13.64% |
| INTRON | 1-100 | POS | 0/0 | 3,697 | 42,501 | 7,389 | 17.39% |
| EXON | 5k-10k | AMBIG | 0/1 | 57 | 20,136 | 7,307 | 36.29% |
| INTRON | 1-100 | NEG | 0/0 | 3,807 | 41,740 | 6,704 | 16.06% |

## Artifacts

- Exact mask confusion: `results/vcap_gdna_false_rna_by_region_kappa_units_fix_2026-05-19/regional_confusion_by_mask.tsv`
- Dominant type confusion: `results/vcap_gdna_false_rna_by_region_kappa_units_fix_2026-05-19/dominant_region_type_confusion.tsv`
- gDNA false RNA by mask: `results/vcap_gdna_false_rna_by_region_kappa_units_fix_2026-05-19/gdna_false_rna_by_mask.tsv`
- Top dominant regions: `results/vcap_gdna_false_rna_by_region_kappa_units_fix_2026-05-19/dominant_region_gdna_false_rna.tsv`
- Dominant region bins: `results/vcap_gdna_false_rna_by_region_kappa_units_fix_2026-05-19/dominant_region_bins.tsv`
- Summary JSON: `results/vcap_gdna_false_rna_by_region_kappa_units_fix_2026-05-19/summary.json`
