# VCaP gDNA False-RNA Hotspot Analysis - 2026-05-16

BAM: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/with_mm/annotated.bam`

Truth source is derived from query-name flowcell ID. This analysis focuses on gDNA-source fragments (`H7MFFDSXY`) that Rigel assigned to an RNA pool (`mRNA` or `nRNA`). Counting rule is one primary read1 record per fragment.

Fragments counted: 32,218,601

gDNA-source fragments: 18,228,677

gDNA -> RNA false positives: 2,912,432 (15.98% of gDNA-source fragments)

False-positive split: 1,340,635 mRNA (46.03% of false RNA) and 1,571,797 nRNA (53.97% of false RNA).

## What Dominates

The false-positive RNA problem is concentrated in unspliced, genic ambiguous fragments. Across false RNA calls, 96.67% are `ZS=unspliced` and 97.30% are either `ambig_same_strand` or `ambig_opp_strand`.

The top 10 EM loci account for 13.62% of all gDNA -> RNA false positives. The top 10 genomic 10,000 bp windows account for 1.09%. So this is not only a few pathological coordinates; it is a recurrent ambiguous-region failure mode, with some strong hotspots.

## Error Categories

| Pred detail | ZC | ZS | Count | False RNA frac | gDNA source frac |
| --- | --- | --- | --- | --- | --- |
| mrna | ambig_same_strand | unspliced | 1,056,095 | 36.26% | 5.79% |
| nrna | ambig_same_strand | unspliced | 1,009,703 | 34.67% | 5.54% |
| nrna | ambig_opp_strand | unspliced | 446,849 | 15.34% | 2.45% |
| mrna | ambig_opp_strand | unspliced | 248,193 | 8.52% | 1.36% |
| nrna | ambig_same_strand | spliced_implicit | 41,963 | 1.44% | 0.23% |
| nrna | unambig | unspliced | 32,188 | 1.11% | 0.18% |
| nrna | multimapper | unspliced | 15,303 | 0.53% | 0.08% |
| nrna | multimapper | spliced_annot | 10,570 | 0.36% | 0.06% |
| mrna | ambig_same_strand | spliced_annot | 9,159 | 0.31% | 0.05% |
| mrna | multimapper | unspliced | 7,134 | 0.24% | 0.04% |
| nrna | ambig_opp_strand | spliced_implicit | 6,826 | 0.23% | 0.04% |
| mrna | ambig_same_strand | spliced_implicit | 6,556 | 0.23% | 0.04% |
| mrna | multimapper | spliced_annot | 4,872 | 0.17% | 0.03% |
| nrna | ambig_same_strand | spliced_unannot | 4,107 | 0.14% | 0.02% |
| nrna | multimapper | spliced_unannot | 3,594 | 0.12% | 0.02% |
| mrna | unambig | spliced_annot | 3,506 | 0.12% | 0.02% |
| mrna | ambig_same_strand | spliced_unannot | 1,850 | 0.06% | 0.01% |
| mrna | ambig_opp_strand | spliced_implicit | 1,664 | 0.06% | 0.01% |
| mrna | multimapper | spliced_unannot | 1,441 | 0.05% | 0.01% |
| nrna | ambig_opp_strand | spliced_unannot | 541 | 0.02% | 0.00% |
| mrna | ambig_opp_strand | spliced_unannot | 165 | 0.01% | 0.00% |
| nrna | unambig | spliced_unannot | 153 | 0.01% | 0.00% |

## Top 10kb Genomic Windows

| Region | False RNA | False rate | gDNA total | mRNA | nRNA | Top target | Top category |
| --- | --- | --- | --- | --- | --- | --- | --- |
| chr11:65,500,001-65,510,000 | 4,011 | 95.36% | 4,206 | 2,008 | 2,003 | nRNA ENST00000698129.1 | nrna/ambig_opp_strand/unspliced |
| chr1:152,350,001-152,360,000 | 3,796 | 82.34% | 4,610 | 1,740 | 2,056 | mRNA FLG2 ENST00000388718.5 | nrna/ambig_opp_strand/unspliced |
| chr11:64,220,001-64,230,000 | 3,350 | 57.67% | 5,809 | 860 | 2,490 | nRNA chr11:64,218,896-64,223,857(+) | nrna/ambig_same_strand/unspliced |
| chr2:178,560,001-178,570,000 | 3,302 | 83.62% | 3,949 | 276 | 3,026 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chrX:67,540,001-67,550,000 | 3,290 | 49.96% | 6,585 | 3,261 | 29 | mRNA AR ENST00000396043.4 | mrna/ambig_same_strand/unspliced |
| chr2:178,570,001-178,580,000 | 3,082 | 82.92% | 3,717 | 210 | 2,872 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr2:178,530,001-178,540,000 | 3,033 | 82.69% | 3,668 | 365 | 2,668 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr11:64,230,001-64,240,000 | 2,750 | 59.56% | 4,617 | 184 | 2,566 | nRNA chr11:64,229,213-64,234,292(-) | nrna/ambig_opp_strand/unspliced |
| chr2:178,540,001-178,550,000 | 2,615 | 83.07% | 3,148 | 103 | 2,512 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr9:106,920,001-106,930,000 | 2,598 | 54.15% | 4,798 | 2,475 | 123 | mRNA ZNF462 ENST00000277225.10 | mrna/ambig_same_strand/unspliced |
| chr2:178,610,001-178,620,000 | 2,192 | 82.53% | 2,656 | 12 | 2,180 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr9:32,630,001-32,640,000 | 2,097 | 60.41% | 3,471 | 348 | 1,749 | nRNA ENST00000242310.4 | nrna/unambig/unspliced |
| chr2:178,590,001-178,600,000 | 2,068 | 82.19% | 2,516 | 73 | 1,995 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr2:178,550,001-178,560,000 | 2,002 | 83.56% | 2,396 | 71 | 1,931 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr21:33,550,001-33,560,000 | 2,002 | 48.57% | 4,122 | 1,884 | 118 | mRNA SON ENST00000300278.8 | mrna/ambig_same_strand/unspliced |
| chr9:76,700,001-76,710,000 | 1,923 | 86.66% | 2,219 | 1,030 | 893 | nRNA chr9:76,704,847-76,710,661(+) | mrna/ambig_opp_strand/unspliced |
| chr17:28,610,001-28,620,000 | 1,918 | 71.59% | 2,679 | 48 | 1,870 | nRNA chr17:28,611,931-28,617,377(+) | nrna/ambig_opp_strand/unspliced |
| chr11:64,110,001-64,120,000 | 1,914 | 47.12% | 4,062 | 1,414 | 500 | mRNA FLRT1 ENST00000246841.3 | mrna/ambig_opp_strand/unspliced |
| chr2:178,580,001-178,590,000 | 1,849 | 80.78% | 2,289 | 52 | 1,797 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr2:178,600,001-178,610,000 | 1,848 | 80.88% | 2,285 | 24 | 1,824 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr11:66,330,001-66,340,000 | 1,840 | 57.37% | 3,207 | 73 | 1,767 | nRNA chr11:66,334,493-66,339,875(+) | nrna/ambig_opp_strand/unspliced |
| chr11:64,260,001-64,270,000 | 1,834 | 34.04% | 5,387 | 58 | 1,776 | nRNA chr11:64,251,529-64,267,923(+) | nrna/ambig_same_strand/unspliced |
| chr2:178,770,001-178,780,000 | 1,779 | 80.10% | 2,221 | 66 | 1,713 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr2:185,790,001-185,800,000 | 1,752 | 83.39% | 2,101 | 855 | 897 | mRNA FSIP2 ENST00000424728.6 | nrna/ambig_opp_strand/unspliced |
| chr15:89,620,001-89,630,000 | 1,734 | 80.43% | 2,156 | 597 | 1,137 | nRNA chr15:89,608,788-89,628,780(-) | nrna/ambig_opp_strand/unspliced |
| chr2:178,720,001-178,730,000 | 1,719 | 75.03% | 2,291 | 8 | 1,711 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr11:66,850,001-66,860,000 | 1,717 | 48.01% | 3,576 | 665 | 1,052 | nRNA chr11:66,857,882-66,860,475(-) | nrna/ambig_opp_strand/unspliced |
| chr2:178,740,001-178,750,000 | 1,686 | 78.49% | 2,148 | 40 | 1,646 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr13:109,780,001-109,790,000 | 1,685 | 43.53% | 3,871 | 1,685 | 0 | mRNA IRS2 ENST00000375856.5 | mrna/ambig_same_strand/unspliced |
| chr2:178,730,001-178,740,000 | 1,664 | 77.22% | 2,155 | 22 | 1,642 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |

## Local RNA Evidence Stratification

False RNA calls are not confined to windows with abundant RNA-source reads. This matters
because it separates two mechanisms: real local RNA expression can pull ambiguous gDNA
into RNA, but ambiguous gDNA can also self-seed RNA components in windows with little or
no RNA evidence.

| RNA-source fragments in 10kb window | Windows | False RNA | False RNA frac | Weighted false rate |
| --- | --- | --- | --- | --- |
| 0 | 17,432 | 376,841 | 12.94% | 14.31% |
| 1-100 | 28,076 | 880,835 | 30.24% | 15.46% |
| 101-1,000 | 16,550 | 1,063,577 | 36.52% | 18.78% |
| >1,000 | 2,975 | 591,179 | 20.30% | 25.28% |

## Top EM Loci

The region span is the observed fragment span among false-positive fragments in the locus, not the full locus extent. Components such as `chr2:3`, `chr1:3`, and `chr7:3` are chromosome-scale mega-components rather than precise biological loci, so the 10kb windows and target tables above are the sharper coordinates for assay triage.

| Locus | Observed region | False RNA | False frac | mRNA | nRNA | Same | Opp | Mean ZW | ZW >= .90 | Top target |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr2:3 | chr2:3,130,406-241,979,860 | 116,670 | 4.01% | 27,518 | 89,152 | 55,511 | 59,816 | 0.565 | 18.19% | nRNA chr2:178,525,988-178,804,642(-) |
| chr1:3 | chr1:21,214-248,936,807 | 45,786 | 1.57% | 18,082 | 27,704 | 28,636 | 14,086 | 0.608 | 17.05% | nRNA chr1:152,189,332-152,365,557(+) |
| chr11:3 | chr11:123,389-130,295,870 | 36,334 | 1.25% | 11,569 | 24,765 | 19,547 | 16,431 | 0.590 | 20.82% | nRNA ENST00000698129.1 |
| chr3:3 | chr3:1,092,672-198,225,656 | 34,020 | 1.17% | 11,088 | 22,932 | 24,096 | 9,683 | 0.576 | 13.37% | nRNA chr3:38,039,204-38,124,025(+) |
| chr7:3 | chr7:997,344-158,587,769 | 33,767 | 1.16% | 11,761 | 22,006 | 24,608 | 8,639 | 0.528 | 10.95% | nRNA chr7:141,995,878-142,106,747(+) |
| chr5:3 | chr5:5,422,760-181,331,826 | 30,985 | 1.06% | 9,697 | 21,288 | 16,604 | 13,328 | 0.630 | 25.17% | nRNA chr5:141,136,682-141,241,744(-) |
| chr9:3 | chr9:14,524-138,331,545 | 26,883 | 0.92% | 11,829 | 15,054 | 17,529 | 7,886 | 0.626 | 23.18% | mRNA ZNF462 ENST00000277225.10 |
| chr8:3 | chr8:120,387-140,458,478 | 26,057 | 0.89% | 8,375 | 17,682 | 18,616 | 7,257 | 0.577 | 15.93% | nRNA chr8:109,362,460-109,537,207(+) |
| chr12:3 | chr12:15,172-133,204,879 | 24,192 | 0.83% | 5,434 | 18,758 | 18,243 | 5,652 | 0.545 | 8.25% | nRNA chr12:40,393,447-40,570,646(+) |
| chr17:3 | chr17:60,017-83,221,642 | 22,110 | 0.76% | 4,521 | 17,589 | 11,069 | 10,524 | 0.584 | 16.67% | nRNA chr17:28,614,881-28,645,454(-) |
| chr4:3 | chr4:1,900,545-190,195,896 | 21,700 | 0.75% | 9,311 | 12,389 | 13,343 | 8,140 | 0.576 | 18.27% | nRNA chr4:78,057,322-78,544,269(+) |
| chr15:3 | chr15:20,144,097-101,968,218 | 18,722 | 0.64% | 4,716 | 14,006 | 9,461 | 9,101 | 0.600 | 15.44% | nRNA chr15:81,331,216-81,374,213(-) |
| chr6:3 | chr6:95,134-170,745,948 | 18,339 | 0.63% | 6,331 | 12,008 | 13,515 | 4,521 | 0.531 | 10.99% | nRNA chr6:56,457,986-56,642,996(-) |
| chr19:3 | chr19:156,517-58,599,378 | 12,056 | 0.41% | 4,832 | 7,224 | 4,863 | 7,057 | 0.622 | 18.08% | nRNA chr19:16,892,950-17,026,815(-) |
| chr13:5941 | chr13:96,090,628-110,512,909 | 9,849 | 0.34% | 448 | 9,401 | 6,780 | 2,976 | 0.557 | 16.15% | nRNA chr13:110,149,060-110,299,118(-) |
| chr5:1380 | chr5:62,306,164-150,396,827 | 8,871 | 0.30% | 1,045 | 7,826 | 8,108 | 693 | 0.739 | 46.80% | nRNA ENST00000617094.1 |
| chr7:883 | chr7:42,757,557-156,969,973 | 8,467 | 0.29% | 2,804 | 5,663 | 4,023 | 3,855 | 0.550 | 10.72% | nRNA chr7:100,197,473-100,214,387(+) |
| chr16:3 | chr16:1,510,646-90,176,006 | 7,214 | 0.25% | 3,264 | 3,950 | 4,724 | 2,083 | 0.552 | 19.98% | mRNA ZFHX3 ENST00000268489.10 |
| chr10:3 | chr10:5,765,702-133,468,331 | 7,138 | 0.25% | 2,930 | 4,208 | 6,379 | 571 | 0.616 | 22.46% | mRNA ZNF518A ENST00000624776.4 |
| chr17:867 | chr17:29,562,271-30,185,966 | 6,414 | 0.22% | 1,699 | 4,715 | 1,624 | 4,776 | 0.599 | 13.42% | nRNA chr17:29,560,546-29,707,090(+) |
| chr19:2804 | chr19:36,054,880-55,456,587 | 5,474 | 0.19% | 1,530 | 3,944 | 4,622 | 710 | 0.480 | 5.01% | nRNA chr19:42,325,634-42,378,769(+) |
| chrX:3 | chrX:7,049,962-155,568,141 | 5,350 | 0.18% | 2,034 | 3,316 | 3,475 | 1,752 | 0.592 | 18.11% | mRNA FRMPD4 ENST00000657982.1 |
| chr14:3 | chr14:16,373,182-81,987,869 | 5,186 | 0.18% | 2,452 | 2,734 | 2,729 | 2,145 | 0.625 | 26.61% | mRNA ADAM20 ENST00000256389.5 |
| chr11:3225 | chr11:5,248,270-64,321,809 | 5,179 | 0.18% | 1,382 | 3,797 | 2,963 | 2,194 | 0.728 | 40.51% | nRNA chr11:64,318,120-64,321,811(+) |
| chr3:1097 | chr3:49,940,497-190,656,877 | 4,835 | 0.17% | 850 | 3,985 | 4,260 | 552 | 0.535 | 12.18% | nRNA chr3:52,495,337-52,524,495(+) |
| chr17:247 | chr17:15,995,681-81,720,975 | 4,654 | 0.16% | 1,108 | 3,546 | 3,773 | 832 | 0.531 | 6.34% | nRNA chr17:81,395,456-81,466,331(+) |
| chr11:3926 | chr11:67,351,579-85,758,027 | 4,574 | 0.16% | 1,243 | 3,331 | 2,511 | 2,043 | 0.506 | 2.58% | nRNA chr11:67,400,854-67,421,183(-) |
| chr18:3 | chr18:1,272,287-65,880,930 | 4,349 | 0.15% | 1,829 | 2,520 | 2,095 | 2,212 | 0.640 | 27.89% | mRNA ALPK2 ENST00000361673.4 |
| chr17:8543 | chr17:10,300,864-46,371,505 | 4,207 | 0.14% | 133 | 4,074 | 23 | 4,160 | 0.533 | 4.61% | nRNA chr17:10,383,155-10,537,862(+) |
| chrX:21523 | chrX:67,521,831-67,853,372 | 4,203 | 0.14% | 3,973 | 230 | 4,166 | 1 | 0.358 | 9.59% | mRNA AR ENST00000396043.4 |

## Top Assigned RNA Targets

| Target | Observed region | False RNA | False frac | mRNA | nRNA | Same | Opp | Mean ZW |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| nRNA chr2:178,525,988-178,804,642(-) | chr2:178,526,713-178,804,640 | 17,642 | 0.61% | 0 | 17,642 | 533 | 17,082 | 0.809 |
| nRNA chr2:178,523,417-178,776,534(+) | chr2:178,526,721-178,776,511 | 9,890 | 0.34% | 0 | 9,890 | 0 | 9,870 | 0.571 |
| nRNA chr11:64,251,529-64,267,923(+) | chr11:64,251,522-64,328,765 | 3,300 | 0.11% | 0 | 3,300 | 3,291 | 0 | 0.810 |
| nRNA chr2:151,485,335-151,734,487(-) | chr2:151,485,477-151,733,424 | 2,477 | 0.09% | 0 | 2,477 | 2,269 | 203 | 0.648 |
| nRNA chr13:110,149,060-110,299,118(-) | chr13:110,149,061-110,298,934 | 2,476 | 0.09% | 0 | 2,476 | 2,467 | 0 | 0.600 |
| nRNA chr2:178,528,764-178,620,148(+) | chr2:32,916,265-178,620,131 | 2,451 | 0.08% | 0 | 2,451 | 0 | 2,448 | 0.200 |
| nRNA chr8:109,362,460-109,537,207(+) | chr8:109,362,486-109,530,370 | 2,371 | 0.08% | 0 | 2,371 | 2,357 | 0 | 0.762 |
| nRNA chr17:29,560,546-29,707,090(+) | chr17:29,562,307-29,703,260 | 2,206 | 0.08% | 0 | 2,206 | 271 | 1,929 | 0.599 |
| mRNA ZFHX4 ENST00000518282.5 | chr8:76,704,046-92,818,668 | 2,190 | 0.08% | 2,190 | 0 | 2,189 | 1 | 0.779 |
| nRNA chr12:57,128,482-57,213,361(+) | chr12:57,128,772-57,212,654 | 2,155 | 0.07% | 0 | 2,155 | 2,064 | 89 | 0.627 |
| nRNA chr1:152,189,332-152,365,557(+) | chr1:152,194,565-152,358,977 | 2,151 | 0.07% | 0 | 2,151 | 12 | 1,662 | 0.670 |
| nRNA chr5:141,136,682-141,241,744(-) | chr5:141,137,134-141,241,485 | 2,146 | 0.07% | 0 | 2,146 | 5 | 2,121 | 0.741 |
| mRNA AHNAK ENST00000378024.9 | chr11:62,516,540-62,534,074 | 2,143 | 0.07% | 2,143 | 0 | 2,130 | 0 | 0.844 |
| mRNA FSIP2 ENST00000424728.6 | chr2:185,788,656-185,814,037 | 2,052 | 0.07% | 2,052 | 0 | 745 | 1,303 | 0.863 |
| nRNA ENST00000698129.1 | chr11:65,499,302-65,506,876 | 1,989 | 0.07% | 0 | 1,989 | 0 | 1,979 | 0.902 |
| mRNA ZNF462 ENST00000277225.10 | chr9:29,052,401-107,011,466 | 1,959 | 0.07% | 1,959 | 0 | 1,936 | 21 | 0.738 |
| mRNA PCLO ENST00000423517.6 | chr7:82,821,985-83,162,894 | 1,959 | 0.07% | 1,959 | 0 | 1,949 | 0 | 0.793 |
| nRNA chr7:141,995,878-142,106,747(+) | chr7:142,005,153-142,106,213 | 1,956 | 0.07% | 0 | 1,956 | 1,930 | 0 | 0.759 |
| nRNA chr11:64,823,051-64,844,653(-) | chr11:64,687,206-64,844,653 | 1,949 | 0.07% | 0 | 1,949 | 1,947 | 0 | 0.800 |
| mRNA MUC16 ENST00000397910.8 | chr19:8,934,899-8,981,316 | 1,852 | 0.06% | 1,852 | 0 | 1,844 | 0 | 0.562 |
| nRNA chr21:8,380,664-8,450,668(+) | mixed:101,576-45,024,528 | 1,815 | 0.06% | 0 | 1,815 | 143 | 0 | 0.827 |
| nRNA chr12:40,393,447-40,570,646(+) | chr12:40,393,447-40,570,627 | 1,788 | 0.06% | 0 | 1,788 | 1,235 | 548 | 0.647 |
| mRNA AR ENST00000396043.4 | chrX:67,521,831-67,724,284 | 1,761 | 0.06% | 1,761 | 0 | 1,747 | 1 | 0.472 |
| mRNA FLG2 ENST00000388718.5 | chr1:152,350,091-152,357,645 | 1,740 | 0.06% | 1,740 | 0 | 0 | 1,725 | 0.791 |
| mRNA IRS2 ENST00000375856.5 | chr13:109,755,789-109,786,266 | 1,736 | 0.06% | 1,736 | 0 | 1,510 | 226 | 0.874 |
| mRNA BSN ENST00000296452.5 | chr3:49,624,972-49,671,352 | 1,679 | 0.06% | 1,679 | 0 | 1,675 | 0 | 0.849 |
| nRNA ENST00000242310.4 | chr9:32,629,863-32,635,660 | 1,643 | 0.06% | 0 | 1,643 | 0 | 644 | 0.891 |
| mRNA FAT1 ENST00000441802.7 | chr4:83,052,074-186,709,839 | 1,583 | 0.05% | 1,583 | 0 | 1,579 | 1 | 0.848 |
| nRNA chr9:110,365,247-110,579,741(-) | chr9:110,365,919-110,579,718 | 1,482 | 0.05% | 0 | 1,482 | 1,480 | 0 | 0.589 |
| nRNA chr11:64,229,213-64,234,292(-) | chr11:64,229,217-64,234,291 | 1,477 | 0.05% | 0 | 1,477 | 0 | 1,470 | 0.829 |

## Representative False-Positive Fragments From Top Loci

| Locus | QNAME | Region | Detail | CIGAR | TLEN | NH | NM | ZW | ZC | ZS | Target |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr2:3 | A00839:75:H7MFFDSXY:3:2324:24279:31955 | chr2:202,765,716-202,765,969 | mrna | 150M | -254 | 1 | 0 | 0.9925652146339417 | ambig_same_strand | unspliced | mRNA FAM117B ENST00000392238.3 |
| chr2:3 | A00839:75:H7MFFDSXY:3:2529:25898:16501 | chr2:151,854,950-151,855,284 | nrna | 150M | 335 | 1 | 0 | 0.07244270294904709 | ambig_same_strand | unspliced | nRNA chr2:151,836,812-151,861,729(-) |
| chr2:3 | A00839:75:H7MFFDSXY:3:1527:21423:23688 | chr2:173,018,541-173,018,752 | nrna | 150M | -212 | 1 | 0 | 0.39551421999931335 | ambig_same_strand | unspliced | nRNA chr2:172,928,204-173,052,893(+) |
| chr2:3 | A00839:75:H7MFFDSXY:3:2602:24777:3850 | chr2:127,261,277-127,261,578 | nrna | 150M | 302 | 1 | 0 | 0.5892437696456909 | ambig_same_strand | unspliced | nRNA chr2:127,257,289-127,262,882(-) |
| chr1:3 | A00839:75:H7MFFDSXY:3:1376:18304:15781 | chr1:109,750,676-109,750,856 | nrna | 150M | 181 | 1 | 0 | 0.7685539722442627 | ambig_opp_strand | unspliced | nRNA chr1:109,750,355-109,751,546(-) |
| chr1:3 | A00839:75:H7MFFDSXY:3:2568:5376:18192 | chr1:33,662,718-33,662,961 | nrna | 150M | 244 | 1 | 0 | 0.181776762008667 | ambig_same_strand | unspliced | nRNA chr1:33,591,806-33,709,504(-) |
| chr1:3 | A00839:75:H7MFFDSXY:3:2130:1289:21762 | chr1:152,353,387-152,353,563 | nrna | 150M | -177 | 1 | 1 | 0.7243297100067139 | ambig_opp_strand | unspliced | nRNA chr1:152,335,378-152,364,080(+) |
| chr1:3 | A00839:75:H7MFFDSXY:3:2674:3161:10441 | chr1:152,351,717-152,351,998 | nrna | 150M | -282 | 1 | 0 | 0.7151299715042114 | ambig_opp_strand | unspliced | nRNA chr1:152,335,378-152,364,080(+) |
| chr11:3 | A00839:75:H7MFFDSXY:3:2218:30761:29277 | chr11:57,737,560-57,737,844 | nrna | 150M | -285 | 1 | 0 | 0.6202215552330017 | ambig_same_strand | unspliced | nRNA chr11:57,712,579-57,740,929(+) |
| chr11:3 | A00839:75:H7MFFDSXY:3:2224:27778:12054 | chr11:65,223,531-65,223,708 | mrna | 150M | -178 | 1 | 0 | 0.5090194344520569 | ambig_same_strand | unspliced | mRNA ENSG00000293475 ENST00000534471.1 |
| chr11:3 | A00839:75:H7MFFDSXY:3:1403:18936:15812 | chr11:65,320,743-65,321,033 | mrna | 150M | 291 | 1 | 0 | 0.1666097491979599 | ambig_same_strand | unspliced | mRNA CDC42EP2 ENST00000533419.1 |
| chr11:3 | A00839:75:H7MFFDSXY:3:1234:31620:12477 | chr11:64,219,947-64,220,282 | nrna | 150M | -336 | 1 | 0 | 0.33219262957572937 | ambig_same_strand | unspliced | nRNA chr11:64,219,452-64,220,482(+) |
| chr3:3 | A00839:75:H7MFFDSXY:3:2437:16893:33786 | chr3:47,343,515-47,343,751 | nrna | 150M | -237 | 1 | 0 | 0.7314984202384949 | ambig_same_strand | unspliced | nRNA chr3:47,282,943-47,344,577(+) |
| chr3:3 | A00839:75:H7MFFDSXY:3:1233:16206:23296 | chr3:53,800,253-53,800,507 | nrna | 150M | -255 | 1 | 0 | 0.9122974872589111 | ambig_same_strand | unspliced | nRNA chr3:53,799,890-53,813,712(+) |
| chr3:3 | A00839:75:H7MFFDSXY:3:2552:23764:25801 | chr3:47,891,475-47,891,772 | mrna | 150M | 298 | 1 | 1 | 0.893859326839447 | ambig_same_strand | unspliced | mRNA MAP4 ENST00000497735.1 |
| chr3:3 | A00839:75:H7MFFDSXY:3:2228:28203:20901 | chr3:38,376,887-38,377,109 | nrna | 150M | -223 | 1 | 0 | 0.8408316969871521 | ambig_same_strand | unspliced | nRNA chr3:38,346,807-38,413,289(+) |
| chr7:3 | A00839:75:H7MFFDSXY:3:2571:8169:34898 | chr7:117,297,519-117,297,759 | nrna | 150M | 241 | 1 | 0 | 0.4051361382007599 | ambig_opp_strand | unspliced | nRNA chr7:117,277,530-117,323,114(-) |
| chr7:3 | A00839:75:H7MFFDSXY:3:1309:10086:6308 | chr7:108,232,308-108,232,482 | nrna | 150M | 175 | 1 | 1 | 0.5315518379211426 | ambig_same_strand | unspliced | nRNA chr7:108,197,994-108,243,240(-) |
| chr7:3 | A00839:75:H7MFFDSXY:3:1666:13919:21277 | chr7:90,163,010-90,163,350 | nrna | 150M | -341 | 1 | 0 | 0.7251846194267273 | ambig_opp_strand | unspliced | nRNA chr7:90,161,265-90,164,829(+) |
| chr7:3 | A00839:75:H7MFFDSXY:3:1155:15709:16986 | chr7:141,919,494-141,919,677 | nrna | 150M | -184 | 1 | 0 | 0.9396408796310425 | ambig_same_strand | unspliced | nRNA ENST00000548136.1 |
| chr5:3 | A00839:75:H7MFFDSXY:3:2337:24813:15843 | chr5:76,295,586-76,295,866 | nrna | 150M | -281 | 1 | 0 | 0.22429996728897095 | ambig_opp_strand | unspliced | nRNA chr5:76,291,764-76,296,122(+) |
| chr5:3 | A00839:75:H7MFFDSXY:3:1134:16306:3114 | chr5:141,200,762-141,201,109 | nrna | 150M | 348 | 1 | 0 | 0.13361455500125885 | ambig_opp_strand | unspliced | nRNA chr5:141,183,400-141,201,396(-) |
| chr5:3 | A00839:75:H7MFFDSXY:3:1163:14262:25441 | chr5:141,183,091-141,183,317 | nrna | 150M | 227 | 1 | 0 | 0.9357520341873169 | ambig_opp_strand | unspliced | nRNA chr5:141,136,682-141,241,744(-) |
| chr5:3 | A00839:75:H7MFFDSXY:3:1512:7401:2409 | chr5:141,235,450-141,235,769 | nrna | 150M | -320 | 1 | 0 | 0.6592143177986145 | ambig_opp_strand | unspliced | nRNA ENST00000524813.2 |
| chr9:3 | A00839:75:H7MFFDSXY:3:1143:16740:3208 | chr9:120,471,757-120,471,980 | nrna | 150M | 224 | 1 | 0 | 0.7240217328071594 | ambig_same_strand | unspliced | nRNA chr9:120,388,874-120,475,733(-) |
| chr9:3 | A00839:75:H7MFFDSXY:3:1647:25554:35289 | chr9:32,633,742-32,634,092 | nrna | 150M | 351 | 1 | 0 | 0.9324710369110107 | ambig_opp_strand | unspliced | nRNA ENST00000242310.4 |
| chr9:3 | A00839:75:H7MFFDSXY:3:2504:2248:3568 | chr9:92,615,947-92,616,184 | mrna | 150M | -238 | 1 | 0 | 0.9618238210678101 | ambig_opp_strand | unspliced | mRNA CENPP ENST00000375587.8 |
| chr9:3 | A00839:75:H7MFFDSXY:3:1347:8467:36385 | chr9:123,073,493-123,073,802 | nrna | 150M | -310 | 1 | 0 | 0.27598318457603455 | ambig_same_strand | unspliced | nRNA chr9:122,941,019-123,104,866(+) |
| chr8:3 | A00839:75:H7MFFDSXY:3:1517:4110:30279 | chr8:99,891,808-99,892,048 | nrna | 150M | 241 | 1 | 1 | 0.4043673276901245 | ambig_same_strand | unspliced | nRNA chr8:99,874,995-99,893,707(-) |
| chr8:3 | A00839:75:H7MFFDSXY:3:1343:12310:18959 | chr8:72,937,822-72,937,992 | mrna | 150M | -171 | 1 | 0 | 0.990434467792511 | ambig_same_strand | unspliced | mRNA KCNB2 ENST00000523207.2 |

## High-Confidence False-Positive Examples

These examples have the largest `ZW` winner weights among gDNA-source fragments called RNA.

| QNAME | Region | Detail | CIGAR | TLEN | NH | NM | ZW | ZC | ZS | Target |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A00839:75:H7MFFDSXY:3:1145:2320:1344 | chr17:47,819,951-47,820,370 | mrna | 150M | -420 | 1 | 0 | 1.0 | ambig_same_strand | spliced_annot | mRNA OSBPL7 ENST00000585051.1 |
| A00839:75:H7MFFDSXY:3:2130:5999:25567 | chr3:184,352,479-184,353,037 | mrna | 145M222N5M | -559 | 1 | 1 | 1.0 | ambig_same_strand | spliced_annot | mRNA CLCN2 ENST00000491162.1 |
| A00839:75:H7MFFDSXY:3:1108:8576:9204 | chr1:154,601,859-154,602,314 | mrna | 150M | -456 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA ADAR ENST00000680270.2 |
| A00839:75:H7MFFDSXY:3:2106:27145:12023 | chrX:109,375,769-109,376,246 | nrna | 150M | -478 | 1 | 0 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chrX:109,372,905-109,482,086(-) |
| A00839:75:H7MFFDSXY:3:1621:15338:32878 | chr9:17,464,264-17,465,983 | nrna | 5S142M1357N3M | -1720 | 2 | 0 | 1.0 | multimapper | spliced_annot | nRNA chr9:17,135,039-17,503,923(+) |
| A00839:75:H7MFFDSXY:3:2262:14733:20275 | chr8:38,457,354-38,457,581 | mrna | 150M | 228 | 13 | 0 | 1.0 | multimapper | unspliced | mRNA FGFR1 ENST00000649678.1 |
| A00839:75:H7MFFDSXY:3:2214:4761:31031 | chr10:102,869,391-102,869,697 | nrna | 119M36N31M | 307 | 1 | 0 | 1.0 | ambig_same_strand | spliced_unannot | nRNA chr10:102,854,258-102,901,899(+) |
| A00839:75:H7MFFDSXY:3:1669:18629:4225 | chr16:20,372,722-20,373,182 | nrna | 150M | -461 | 1 | 0 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chr16:20,359,174-20,404,737(-) |
| A00839:75:H7MFFDSXY:3:2627:2040:11600 | chr11:65,617,181-65,617,638 | nrna | 150M | 458 | 1 | 0 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chr11:65,615,775-65,637,439(+) |
| A00839:75:H7MFFDSXY:3:1561:32470:6151 | chr2:85,306,158-85,306,633 | nrna | 150M | 476 | 1 | 0 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chr2:85,133,391-85,310,387(+) |
| A00839:75:H7MFFDSXY:3:2358:28456:2519 | chr3:52,451,133-52,451,579 | nrna | 150M | -447 | 1 | 0 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chr3:52,451,099-52,454,041(-) |
| A00839:75:H7MFFDSXY:3:2571:21260:24064 | chr19:2,345,068-2,345,502 | nrna | 150M | 435 | 40 | 1 | 1.0 | multimapper | unspliced | nRNA chr19:24,033,400-24,129,968(+) |
| A00839:75:H7MFFDSXY:3:2418:20880:1454 | chr9:134,317,379-134,317,685 | nrna | 136M26N14M | 307 | 1 | 0 | 1.0 | ambig_same_strand | spliced_unannot | nRNA chr9:134,317,097-134,406,394(+) |
| A00839:75:H7MFFDSXY:3:1627:12292:24940 | chr7:128,837,346-128,837,784 | nrna | 150M | -439 | 1 | 0 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chr7:128,830,405-128,859,274(+) |
| A00839:75:H7MFFDSXY:3:1405:4643:18239 | chrX:41,342,042-41,342,822 | mrna | 150M | -781 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA DDX3X ENST00000615742.4 |
| A00839:75:H7MFFDSXY:3:1550:11704:5102 | chr20:63,542,163-63,542,607 | nrna | 150M | 445 | 1 | 1 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chr20:63,538,488-63,547,749(-) |
| A00839:75:H7MFFDSXY:3:2475:1199:2879 | chr1:156,585,632-156,586,213 | nrna | 150M | 582 | 1 | 0 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chr1:156,579,722-156,587,719(+) |
| A00839:75:H7MFFDSXY:3:2574:6876:15342 | chr14:91,303,617-91,305,767 | nrna | 150M | 2151 | 2 | 0 | 1.0 | multimapper | spliced_annot | nRNA chr14:91,271,322-91,417,820(-) |
| A00839:75:H7MFFDSXY:3:1635:8730:19147 | chr17:8,836,437-8,836,794 | mrna | 146M74N4M | -358 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA PIK3R6 ENST00000611951.4 |
| A00839:75:H7MFFDSXY:3:1322:26232:31548 | chr5:146,878,015-146,878,593 | mrna | 147M395N3M | -579 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA PPP2R2B ENST00000394411.9 |
| A00839:75:H7MFFDSXY:3:1245:16902:7811 | chr1:43,364,150-43,364,617 | mrna | 147M154N3M | -468 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA ELOVL1 ENST00000464204.5 |
| A00839:75:H7MFFDSXY:3:1658:8043:31046 | chr11:126,205,396-126,205,744 | mrna | 144M168N6M | 349 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA RPUSD4 ENST00000533628.5 |
| A00839:75:H7MFFDSXY:3:2604:17996:32503 | chr17:1,745,211-1,745,630 | mrna | 150M | -420 | 1 | 0 | 1.0 | ambig_same_strand | spliced_annot | mRNA SERPINF2 ENST00000450523.6 |
| A00839:75:H7MFFDSXY:3:2418:30572:12414 | chr2:85,133,768-85,134,234 | nrna | 150M | -467 | 1 | 0 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chr2:85,133,391-85,310,387(+) |
| A00839:75:H7MFFDSXY:3:1265:1542:19539 | chr12:120,534,243-120,534,805 | nrna | 150M | 563 | 1 | 0 | 1.0 | ambig_opp_strand | spliced_implicit | nRNA chr12:120,533,479-120,534,836(+) |
| A00839:75:H7MFFDSXY:3:2271:22724:9815 | chr9:35,710,735-35,711,214 | nrna | 150M | 480 | 1 | 0 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chr9:35,696,947-35,732,195(-) |
| A00839:75:H7MFFDSXY:3:1355:5665:1971 | chr22:49,187,997-49,188,386 | nrna | 132M | -390 | 4 | 1 | 1.0 | multimapper | spliced_unannot | nRNA chr22:49,166,420-49,188,062(-) |
| A00839:75:H7MFFDSXY:3:2375:22119:36949 | chr1:47,035,751-47,036,196 | nrna | 150M | 446 | 1 | 0 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chr1:47,023,668-47,050,751(+) |
| A00839:75:H7MFFDSXY:3:1260:29441:17221 | chr12:8,917,724-8,921,043 | mrna | 91M1067N59M | -3320 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA PHC1 ENST00000541181.1 |
| A00839:75:H7MFFDSXY:3:2610:13476:33755 | chr11:64,827,252-64,827,548 | mrna | 147M128N3M | 297 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA CDC42BPG ENST00000342711.6 |

## Interpretation

1. The dominant false positives are not spliced RNA-like evidence. They are ordinary contiguous gDNA alignments (`ZS=unspliced`, usually `NH=1`) that overlap expressed annotated transcript or nRNA spans.
2. Ambiguous same-strand and opposite-strand labels mean the fragment has RNA candidates but no decisive splice-junction or transcript-unique evidence. In these cases, the EM can let RNA abundance in a large component overcome the generic gDNA candidate.
3. Local expression is only part of the failure. Many false-positive windows have little or no RNA-source support, meaning gDNA fragments themselves can self-seed an RNA component when the model has no conservative guard for ambiguous unspliced evidence.
4. Many false positives have high `ZW`, so this is not just low-confidence noise. The model is often confident after EM because the RNA component has enough mass to explain genic unspliced fragments.
5. A conservative diagnostic mode should treat ambiguous unspliced genic fragments from gDNA-compatible regions as weak RNA evidence unless they are supported by splice junctions, transcript-unique exonic geometry, or local RNA-only context.

## Recommended Changes

1. Add an RNA-call guard for `ZC in {ambig_same_strand, ambig_opp_strand}` and `ZS=unspliced`: require a stronger RNA-vs-gDNA posterior margin before emitting mRNA/nRNA as the hard winner.
2. Make gDNA a first-class competitor for ambiguous genic unspliced fragments even when locus RNA mass is high; practically, this means a conservative prior or likelihood floor for gDNA in these channels.
3. Split RNA evidence tiers in diagnostic output: splice-supported RNA, transcript-unique exonic RNA, and ambiguous-unspliced RNA. The last tier should not be used as strong positive RNA evidence in assays without additional support.
4. Preserve enough annotation in BAM output to inspect the runner-up gDNA posterior or RNA/gDNA posterior ratio. `ZW` alone reports winner confidence, but not the rejected gDNA probability.
5. Consider a diagnostic operating mode that defaults ambiguous unspliced RNA/gDNA ties to gDNA, trading RNA sensitivity for much lower false-positive RNA.
6. Break or regularize chromosome-scale mega-components for diagnostic calling. Ambiguous unspliced fragments should be judged against local evidence, not allowed to borrow RNA mass from distant transcripts in the same connected component.

## Artifacts

- Error categories: `results/vcap_with_mm_hotspots_2026-05-17/gdna_false_rna_categories.tsv`
- Hotspot windows: `results/vcap_with_mm_hotspots_2026-05-17/gdna_false_rna_windows.tsv`
- Hotspot EM loci: `results/vcap_with_mm_hotspots_2026-05-17/gdna_false_rna_loci.tsv`
- Assigned RNA targets: `results/vcap_with_mm_hotspots_2026-05-17/gdna_false_rna_targets.tsv`
- Representative fragments: `results/vcap_with_mm_hotspots_2026-05-17/sample_false_rna_fragments.tsv`
- High-confidence examples: `results/vcap_with_mm_hotspots_2026-05-17/high_confidence_false_rna_fragments.tsv`
- Summary JSON: `results/vcap_with_mm_hotspots_2026-05-17/summary.json`
