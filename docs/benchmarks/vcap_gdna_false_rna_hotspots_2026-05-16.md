# VCaP gDNA False-RNA Hotspot Analysis - 2026-05-16

BAM: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam`

Truth source is derived from query-name flowcell ID. This analysis focuses on gDNA-source fragments (`H7MFFDSXY`) that Rigel assigned to an RNA pool (`mRNA` or `nRNA`). Counting rule is one primary read1 record per fragment.

Fragments counted: 32,218,601

gDNA-source fragments: 18,228,677

gDNA -> RNA false positives: 3,331,662 (18.28% of gDNA-source fragments)

False-positive split: 1,385,117 mRNA (41.57% of false RNA) and 1,946,545 nRNA (58.43% of false RNA).

## What Dominates

The false-positive RNA problem is concentrated in unspliced, genic ambiguous fragments. Across false RNA calls, 99.54% are `ZS=unspliced` and 97.85% are either `ambig_same_strand` or `ambig_opp_strand`.

The top 10 EM loci account for 8.58% of all gDNA -> RNA false positives. The top 10 genomic 10,000 bp windows account for 0.65%. So this is not only a few pathological coordinates; it is a recurrent ambiguous-region failure mode, with some strong hotspots.

## Error Categories

| Pred detail | ZC | ZS | Count | False RNA frac | gDNA source frac |
| --- | --- | --- | --- | --- | --- |
| nrna | ambig_same_strand | unspliced | 1,299,008 | 38.99% | 7.13% |
| mrna | ambig_same_strand | unspliced | 1,107,859 | 33.25% | 6.08% |
| nrna | ambig_opp_strand | unspliced | 582,324 | 17.48% | 3.19% |
| mrna | ambig_opp_strand | unspliced | 262,978 | 7.89% | 1.44% |
| nrna | unambig | unspliced | 35,519 | 1.07% | 0.19% |
| nrna | multimapper | unspliced | 20,898 | 0.63% | 0.11% |
| mrna | multimapper | unspliced | 7,705 | 0.23% | 0.04% |
| nrna | ambig_same_strand | spliced_unannot | 3,773 | 0.11% | 0.02% |
| nrna | multimapper | spliced_unannot | 2,694 | 0.08% | 0.01% |
| mrna | ambig_same_strand | spliced_annot | 2,095 | 0.06% | 0.01% |
| nrna | multimapper | spliced_annot | 1,835 | 0.06% | 0.01% |
| mrna | ambig_same_strand | spliced_unannot | 1,712 | 0.05% | 0.01% |
| mrna | multimapper | spliced_unannot | 1,200 | 0.04% | 0.01% |
| mrna | multimapper | spliced_annot | 863 | 0.03% | 0.00% |
| mrna | unambig | spliced_annot | 616 | 0.02% | 0.00% |
| nrna | ambig_opp_strand | spliced_unannot | 352 | 0.01% | 0.00% |
| nrna | unambig | spliced_unannot | 142 | 0.00% | 0.00% |
| mrna | ambig_opp_strand | spliced_unannot | 89 | 0.00% | 0.00% |

## Top 10kb Genomic Windows

| Region | False RNA | False rate | gDNA total | mRNA | nRNA | Top target | Top category |
| --- | --- | --- | --- | --- | --- | --- | --- |
| chr1:152,350,001-152,360,000 | 3,253 | 70.56% | 4,610 | 1,587 | 1,666 | mRNA FLG2 ENST00000388718.5 | nrna/ambig_opp_strand/unspliced |
| chrX:67,540,001-67,550,000 | 2,885 | 43.81% | 6,585 | 2,881 | 4 | mRNA AR ENST00000396043.4 | mrna/ambig_same_strand/unspliced |
| chr9:32,630,001-32,640,000 | 2,245 | 64.68% | 3,471 | 374 | 1,871 | nRNA ENST00000242310.4 | nrna/unambig/unspliced |
| chr9:106,920,001-106,930,000 | 2,101 | 43.79% | 4,798 | 2,040 | 61 | mRNA ZNF462 ENST00000277225.10 | mrna/ambig_same_strand/unspliced |
| chr9:76,700,001-76,710,000 | 2,039 | 91.89% | 2,219 | 1,034 | 1,005 | nRNA chr9:76,704,847-76,710,661(+) | mrna/ambig_opp_strand/unspliced |
| chr11:64,110,001-64,120,000 | 1,978 | 48.70% | 4,062 | 1,469 | 509 | mRNA FLRT1 ENST00000246841.3 | mrna/ambig_opp_strand/unspliced |
| chr15:89,620,001-89,630,000 | 1,922 | 89.15% | 2,156 | 635 | 1,287 | nRNA chr15:89,608,788-89,628,780(-) | nrna/ambig_opp_strand/unspliced |
| chr21:33,550,001-33,560,000 | 1,902 | 46.14% | 4,122 | 1,764 | 138 | mRNA SON ENST00000300278.8 | mrna/ambig_same_strand/unspliced |
| chr11:66,850,001-66,860,000 | 1,750 | 48.94% | 3,576 | 631 | 1,119 | nRNA chr11:66,857,882-66,860,475(-) | nrna/ambig_opp_strand/unspliced |
| chr11:65,500,001-65,510,000 | 1,721 | 40.92% | 4,206 | 1,720 | 1 | mRNA MALAT1 ENST00000618132.2 | mrna/ambig_opp_strand/unspliced |
| chr17:29,570,001-29,580,000 | 1,685 | 57.51% | 2,930 | 134 | 1,551 | nRNA chr17:29,560,546-29,707,090(+) | nrna/ambig_opp_strand/unspliced |
| chr13:109,780,001-109,790,000 | 1,666 | 43.04% | 3,871 | 1,666 | 0 | mRNA IRS2 ENST00000375856.5 | mrna/ambig_same_strand/unspliced |
| chr4:87,610,001-87,620,000 | 1,637 | 72.69% | 2,252 | 751 | 886 | mRNA DSPP ENST00000651931.1 | nrna/ambig_opp_strand/unspliced |
| chr19:56,810,001-56,820,000 | 1,620 | 69.59% | 2,328 | 831 | 789 | nRNA chr19:56,765,318-56,823,395(+) | mrna/ambig_opp_strand/unspliced |
| chr2:185,790,001-185,800,000 | 1,553 | 73.92% | 2,101 | 790 | 763 | mRNA FSIP2 ENST00000424728.6 | mrna/ambig_opp_strand/unspliced |
| chr5:141,170,001-141,180,000 | 1,539 | 87.54% | 1,758 | 0 | 1,539 | nRNA chr5:141,136,682-141,241,744(-) | nrna/ambig_opp_strand/unspliced |
| chr11:67,250,001-67,260,000 | 1,464 | 55.45% | 2,640 | 518 | 946 | nRNA chr11:67,252,335-67,261,545(-) | nrna/ambig_opp_strand/unspliced |
| chr13:23,330,001-23,340,000 | 1,418 | 41.55% | 3,413 | 1,418 | 0 | mRNA SACS ENST00000402364.1 | mrna/ambig_same_strand/unspliced |
| chr8:105,800,001-105,810,000 | 1,397 | 82.66% | 1,690 | 740 | 657 | nRNA chr8:105,785,054-105,911,869(-) | mrna/ambig_opp_strand/unspliced |
| chr17:29,630,001-29,640,000 | 1,383 | 62.69% | 2,206 | 696 | 687 | mRNA SSH2 ENST00000649863.1 | mrna/ambig_opp_strand/unspliced |
| chr8:76,850,001-76,860,000 | 1,376 | 47.12% | 2,920 | 1,342 | 34 | mRNA ZFHX4 ENST00000518282.5 | mrna/ambig_same_strand/unspliced |
| chr17:29,610,001-29,620,000 | 1,362 | 49.85% | 2,732 | 100 | 1,262 | nRNA chr17:29,560,546-29,707,090(+) | nrna/ambig_opp_strand/unspliced |
| chr15:99,130,001-99,140,000 | 1,359 | 62.68% | 2,168 | 922 | 437 | mRNA SYNM ENST00000336292.11 | mrna/ambig_same_strand/unspliced |
| chr8:54,620,001-54,630,000 | 1,314 | 43.17% | 3,044 | 1,132 | 182 | mRNA RP1 ENST00000220676.2 | mrna/ambig_same_strand/unspliced |
| chr11:64,280,001-64,290,000 | 1,311 | 29.40% | 4,459 | 544 | 767 | nRNA chr11:64,286,419-64,289,494(+) | nrna/ambig_same_strand/unspliced |
| chr11:65,860,001-65,870,000 | 1,267 | 26.74% | 4,739 | 133 | 1,134 | nRNA chr11:65,860,061-65,866,440(+) | nrna/ambig_opp_strand/unspliced |
| chr11:64,220,001-64,230,000 | 1,238 | 21.31% | 5,809 | 409 | 829 | nRNA chr11:64,223,798-64,226,158(-) | nrna/ambig_same_strand/unspliced |
| chr11:66,280,001-66,290,000 | 1,237 | 43.82% | 2,823 | 65 | 1,172 | nRNA chr11:66,278,461-66,285,301(+) | nrna/ambig_same_strand/unspliced |
| chr3:49,650,001-49,660,000 | 1,232 | 37.29% | 3,304 | 1,232 | 0 | mRNA BSN ENST00000296452.5 | mrna/ambig_same_strand/unspliced |
| chr5:140,640,001-140,650,000 | 1,225 | 79.39% | 1,543 | 67 | 1,158 | nRNA chr5:140,638,739-140,647,314(-) | nrna/ambig_opp_strand/unspliced |

## Local RNA Evidence Stratification

False RNA calls are not confined to windows with abundant RNA-source reads. This matters
because it separates two mechanisms: real local RNA expression can pull ambiguous gDNA
into RNA, but ambiguous gDNA can also self-seed RNA components in windows with little or
no RNA evidence.

| RNA-source fragments in 10kb window | Windows | False RNA | False RNA frac | Weighted false rate |
| --- | --- | --- | --- | --- |
| 0 | 26,136 | 491,471 | 14.75% | 17.03% |
| 1-100 | 32,207 | 1,084,447 | 32.55% | 18.21% |
| 101-1,000 | 16,928 | 1,186,455 | 35.61% | 20.81% |
| >1,000 | 2,972 | 569,289 | 17.09% | 24.36% |

## Top EM Loci

The region span is the observed fragment span among false-positive fragments in the locus, not the full locus extent. Components such as `chr2:3`, `chr1:3`, and `chr7:3` are chromosome-scale mega-components rather than precise biological loci, so the 10kb windows and target tables above are the sharper coordinates for assay triage.

| Locus | Observed region | False RNA | False frac | mRNA | nRNA | Same | Opp | Mean ZW | ZW >= .90 | Top target |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr2:3 | chr2:3,188,987-231,371,858 | 37,204 | 1.12% | 9,799 | 27,405 | 26,076 | 10,871 | 0.607 | 14.61% | nRNA chr2:151,485,335-151,734,487(-) |
| chr1:3 | chr1:15,551-248,936,807 | 34,842 | 1.05% | 9,913 | 24,929 | 26,684 | 6,185 | 0.605 | 18.15% | nRNA chr1:237,471,173-237,833,988(+) |
| chr7:3 | chr7:997,344-158,587,748 | 34,397 | 1.03% | 8,733 | 25,664 | 24,612 | 9,433 | 0.560 | 11.34% | nRNA chr7:141,995,878-142,106,747(+) |
| chr3:3 | chr3:3,799,531-198,225,656 | 32,759 | 0.98% | 8,602 | 24,157 | 22,009 | 10,630 | 0.608 | 19.70% | nRNA chr3:38,039,204-38,124,025(+) |
| chr4:3 | chr4:9,130,455-190,195,896 | 27,694 | 0.83% | 8,586 | 19,108 | 16,121 | 11,391 | 0.610 | 15.32% | nRNA chr4:78,057,322-78,544,269(+) |
| chr5:3 | chr5:2,087,416-181,331,826 | 26,229 | 0.79% | 5,679 | 20,550 | 11,762 | 14,076 | 0.692 | 32.66% | nRNA chr5:141,136,682-141,241,744(-) |
| chr9:3 | chr9:14,524-138,331,545 | 24,983 | 0.75% | 7,918 | 17,065 | 15,050 | 8,399 | 0.665 | 23.76% | nRNA chr9:110,365,247-110,579,741(-) |
| chr8:3 | chr8:166,492-140,458,521 | 24,173 | 0.73% | 5,327 | 18,846 | 18,225 | 5,866 | 0.604 | 20.54% | nRNA chr8:109,362,460-109,537,207(+) |
| chr15:3 | chr15:20,144,097-101,968,218 | 21,764 | 0.65% | 4,285 | 17,479 | 11,073 | 10,567 | 0.620 | 16.46% | nRNA chr15:81,331,216-81,374,213(-) |
| chr6:3 | chr6:95,134-170,745,960 | 21,759 | 0.65% | 6,047 | 15,712 | 16,572 | 5,018 | 0.559 | 13.35% | nRNA chr6:56,457,986-56,642,996(-) |
| chr12:3 | chr12:15,091-133,204,855 | 20,887 | 0.63% | 3,663 | 17,224 | 16,077 | 4,658 | 0.580 | 10.64% | nRNA chr12:26,335,351-26,833,194(-) |
| chr11:3 | chr11:87,377-130,301,460 | 16,388 | 0.49% | 4,456 | 11,932 | 10,714 | 5,492 | 0.603 | 17.01% | nRNA chr11:46,856,716-46,918,550(-) |
| chr17:3 | chr17:60,017-83,221,780 | 12,792 | 0.38% | 3,259 | 9,533 | 6,972 | 5,674 | 0.546 | 10.92% | nRNA chr17:31,258,382-31,336,473(+) |
| chr19:3 | chr19:61,553-58,583,462 | 11,647 | 0.35% | 3,583 | 8,064 | 3,830 | 7,736 | 0.679 | 26.61% | nRNA chr19:16,892,950-17,026,815(-) |
| chr7:917 | chr7:43,926,439-102,808,838 | 9,591 | 0.29% | 2,126 | 7,465 | 3,669 | 5,172 | 0.581 | 10.85% | nRNA chr7:100,197,473-100,214,387(+) |
| chr16:3 | chr16:1,510,646-90,176,006 | 8,055 | 0.24% | 2,510 | 5,545 | 5,011 | 2,521 | 0.575 | 20.29% | mRNA ZFHX3 ENST00000268489.10 |
| chr10:3 | chr10:5,765,702-133,468,331 | 7,693 | 0.23% | 2,257 | 5,436 | 7,373 | 202 | 0.631 | 17.87% | nRNA chr10:99,782,639-99,852,594(+) |
| chr1:1344 | chr1:152,124,317-152,412,487 | 7,445 | 0.22% | 4,455 | 2,990 | 904 | 5,748 | 0.617 | 10.44% | mRNA FLG2 ENST00000388718.5 |
| chr2:12673 | chr2:178,526,893-178,802,535 | 7,438 | 0.22% | 826 | 6,612 | 56 | 7,371 | 0.150 | 0.07% | nRNA chr2:178,525,988-178,804,642(-) |
| chrX:3 | chrX:7,049,962-155,524,630 | 7,137 | 0.21% | 2,197 | 4,940 | 4,714 | 2,344 | 0.601 | 18.00% | nRNA chrX:65,499,117-65,534,757(-) |
| chr17:900 | chr17:29,562,177-30,185,957 | 7,130 | 0.21% | 1,760 | 5,370 | 1,736 | 5,384 | 0.617 | 12.73% | nRNA chr17:29,560,546-29,707,090(+) |
| chr17:8885 | chr17:10,300,771-10,729,383 | 6,544 | 0.20% | 129 | 6,415 | 42 | 6,479 | 0.511 | 1.04% | nRNA chr17:10,383,155-10,537,862(+) |
| chr14:3 | chr14:20,510,243-81,987,869 | 5,833 | 0.18% | 2,333 | 3,500 | 2,813 | 2,973 | 0.617 | 23.69% | mRNA ADAM20 ENST00000256389.5 |
| chr19:2934 | chr19:42,325,719-55,374,928 | 5,292 | 0.16% | 1,286 | 4,006 | 4,292 | 948 | 0.479 | 1.53% | nRNA chr19:42,325,634-42,378,769(+) |
| chr5:17546 | chr5:141,330,638-141,512,199 | 4,983 | 0.15% | 416 | 4,567 | 4,784 | 163 | 0.665 | 5.02% | nRNA ENST00000611950.1 |
| chr20:3 | chr20:1,757,716-64,292,881 | 4,536 | 0.14% | 2,211 | 2,325 | 3,495 | 990 | 0.688 | 36.77% | mRNA TSHZ2 ENST00000603338.2 |
| chr11:3369 | chr11:5,248,270-5,904,282 | 4,514 | 0.14% | 1,450 | 3,064 | 1,641 | 2,850 | 0.697 | 39.96% | nRNA ENST00000380184.2 |
| chr2:9930 | chr2:74,455,093-202,560,197 | 4,280 | 0.13% | 1,782 | 2,498 | 3,083 | 1,106 | 0.652 | 33.83% | mRNA RANBP2 ENST00000283195.11 |
| chr17:9038 | chr17:28,473,392-28,645,234 | 4,135 | 0.12% | 897 | 3,238 | 2,202 | 1,927 | 0.537 | 4.30% | nRNA chr17:28,614,881-28,645,454(-) |
| chr21:13954 | chr21:33,503,934-33,968,616 | 4,104 | 0.12% | 2,360 | 1,744 | 3,099 | 988 | 0.476 | 6.04% | mRNA SON ENST00000300278.8 |

## Top Assigned RNA Targets

| Target | Observed region | False RNA | False frac | mRNA | nRNA | Same | Opp | Mean ZW |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| nRNA chr2:151,485,335-151,734,487(-) | chr2:151,485,401-151,733,434 | 3,097 | 0.09% | 0 | 3,097 | 2,857 | 234 | 0.742 |
| nRNA chr8:109,362,460-109,537,207(+) | chr8:109,362,486-109,530,469 | 2,932 | 0.09% | 0 | 2,932 | 2,922 | 0 | 0.837 |
| nRNA chr5:141,136,682-141,241,744(-) | chr5:141,136,761-141,241,568 | 2,639 | 0.08% | 0 | 2,639 | 13 | 2,595 | 0.795 |
| nRNA chr17:29,560,546-29,707,090(+) | chr17:29,562,307-29,703,260 | 2,575 | 0.08% | 0 | 2,575 | 296 | 2,273 | 0.650 |
| nRNA chr7:141,995,878-142,106,747(+) | chr7:142,005,153-142,106,219 | 2,223 | 0.07% | 0 | 2,223 | 2,200 | 0 | 0.779 |
| mRNA AHNAK ENST00000378024.9 | chr11:62,516,540-62,534,074 | 2,210 | 0.07% | 2,210 | 0 | 2,198 | 0 | 0.861 |
| nRNA chr4:78,057,322-78,544,269(+) | mixed:25,647,226-78,541,173 | 2,157 | 0.06% | 0 | 2,157 | 2,150 | 0 | 0.653 |
| nRNA chr9:110,365,247-110,579,741(-) | chr9:110,365,919-110,579,730 | 2,041 | 0.06% | 0 | 2,041 | 2,041 | 0 | 0.694 |
| mRNA MUC16 ENST00000397910.8 | chr19:8,934,968-8,981,316 | 2,025 | 0.06% | 2,025 | 0 | 2,016 | 0 | 0.563 |
| mRNA PCLO ENST00000423517.6 | chr7:82,821,985-83,162,905 | 2,018 | 0.06% | 2,018 | 0 | 2,010 | 0 | 0.804 |
| mRNA ZFHX4 ENST00000518282.5 | chr8:76,704,045-76,864,778 | 2,014 | 0.06% | 2,014 | 0 | 2,014 | 0 | 0.839 |
| mRNA FSIP2 ENST00000424728.6 | chr2:185,788,659-185,833,291 | 1,902 | 0.06% | 1,902 | 0 | 692 | 1,206 | 0.848 |
| nRNA chr21:8,380,664-8,450,668(+) | mixed:101,586-8,472,277 | 1,868 | 0.06% | 0 | 1,868 | 129 | 0 | 0.844 |
| mRNA BSN ENST00000296452.5 | chr3:49,624,975-49,671,352 | 1,737 | 0.05% | 1,737 | 0 | 1,735 | 0 | 0.869 |
| nRNA ENST00000242310.4 | chr9:32,629,863-32,635,660 | 1,722 | 0.05% | 0 | 1,722 | 0 | 658 | 0.951 |
| mRNA IRS2 ENST00000375856.5 | chr13:109,755,796-109,786,266 | 1,715 | 0.05% | 1,715 | 0 | 1,561 | 154 | 0.911 |
| nRNA chr2:178,525,988-178,804,642(-) | chr2:178,526,893-178,802,535 | 1,676 | 0.05% | 0 | 1,676 | 49 | 1,625 | 0.121 |
| nRNA chr12:123,762,187-123,935,720(+) | chr12:123,762,190-123,935,691 | 1,673 | 0.05% | 0 | 1,673 | 1,452 | 215 | 0.596 |
| mRNA FAT3 ENST00000409404.6 | chr11:92,352,096-92,891,554 | 1,666 | 0.05% | 1,666 | 0 | 1,654 | 12 | 0.828 |
| mRNA ZNF462 ENST00000277225.10 | chr9:106,923,398-107,011,446 | 1,654 | 0.05% | 1,654 | 0 | 1,634 | 17 | 0.788 |
| mRNA FAT1 ENST00000441802.7 | chr4:186,588,211-186,709,840 | 1,620 | 0.05% | 1,620 | 0 | 1,618 | 0 | 0.864 |
| nRNA chr6:56,457,986-56,642,996(-) | chr6:56,458,896-56,642,933 | 1,606 | 0.05% | 0 | 1,606 | 1,605 | 0 | 0.493 |
| mRNA AR ENST00000396043.4 | chrX:67,521,831-67,724,284 | 1,588 | 0.05% | 1,588 | 0 | 1,580 | 0 | 0.518 |
| mRNA FLG2 ENST00000388718.5 | chr1:152,350,091-152,357,645 | 1,587 | 0.05% | 1,587 | 0 | 0 | 1,571 | 0.770 |
| nRNA chr9:114,850,967-115,121,030(+) | chr9:114,862,696-115,106,817 | 1,529 | 0.05% | 0 | 1,529 | 16 | 1,507 | 0.778 |
| mRNA FAT4 ENST00000394329.9 | chr4:125,315,073-125,491,980 | 1,491 | 0.04% | 1,491 | 0 | 1,489 | 0 | 0.811 |
| mRNA FLRT1 ENST00000246841.3 | chr11:64,103,236-64,118,848 | 1,469 | 0.04% | 1,469 | 0 | 0 | 1,469 | 0.861 |
| nRNA chr2:169,127,108-169,362,534(-) | chr2:169,128,526-169,362,517 | 1,446 | 0.04% | 0 | 1,446 | 1,445 | 0 | 0.639 |
| mRNA SON ENST00000300278.8 | chr21:33,549,478-33,560,571 | 1,391 | 0.04% | 1,391 | 0 | 1,331 | 59 | 0.650 |
| mRNA SACS ENST00000402364.1 | chr13:23,329,055-23,375,270 | 1,387 | 0.04% | 1,387 | 0 | 1,387 | 0 | 0.674 |

## Representative False-Positive Fragments From Top Loci

| Locus | QNAME | Region | Detail | CIGAR | TLEN | NH | NM | ZW | ZC | ZS | Target |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr2:3 | A00839:75:H7MFFDSXY:3:2324:24279:31955 | chr2:202,765,716-202,765,969 | mrna | 150M | -254 | 1 | 0 | 0.9963996410369873 | ambig_same_strand | unspliced | mRNA FAM117B ENST00000392238.3 |
| chr2:3 | A00839:75:H7MFFDSXY:3:2529:25898:16501 | chr2:151,854,950-151,855,284 | mrna | 150M | 335 | 1 | 0 | 0.6963960528373718 | ambig_same_strand | unspliced | mRNA CACNB4 ENST00000636834.1 |
| chr2:3 | A00839:75:H7MFFDSXY:3:1527:2211:4100 | chr2:167,245,830-167,246,212 | mrna | 150M | -383 | 1 | 0 | 0.19109076261520386 | ambig_same_strand | unspliced | mRNA XIRP2 ENST00000672277.1 |
| chr2:3 | A00839:75:H7MFFDSXY:3:1152:12581:34710 | chr2:151,696,491-151,696,758 | nrna | 150M | 268 | 1 | 0 | 0.91264808177948 | ambig_same_strand | unspliced | nRNA chr2:151,485,335-151,734,487(-) |
| chr1:3 | A00839:75:H7MFFDSXY:3:2525:31439:29262 | chr1:146,685,556-146,685,814 | nrna | 150M | -259 | 2 | 0 | 0.633751630783081 | multimapper | unspliced | nRNA chr1:146,547,366-146,898,974(+) |
| chr1:3 | A00839:75:H7MFFDSXY:3:2370:2510:25755 | chr1:240,143,811-240,144,209 | nrna | 150M | 399 | 1 | 0 | 0.5069261789321899 | ambig_opp_strand | unspliced | nRNA ENST00000451892.1 |
| chr1:3 | A00839:75:H7MFFDSXY:3:1376:18304:15781 | chr1:109,750,676-109,750,856 | nrna | 150M | 181 | 1 | 0 | 0.748022735118866 | ambig_opp_strand | unspliced | nRNA chr1:109,750,355-109,751,546(-) |
| chr1:3 | A00839:75:H7MFFDSXY:3:2568:5376:18192 | chr1:33,662,718-33,662,961 | nrna | 150M | 244 | 1 | 0 | 0.33602648973464966 | ambig_same_strand | unspliced | nRNA chr1:33,513,997-34,165,842(-) |
| chr7:3 | A00839:75:H7MFFDSXY:3:1421:6072:1047 | chr7:117,385,568-117,385,863 | nrna | 150M | 296 | 1 | 0 | 0.6558648347854614 | ambig_opp_strand | unspliced | nRNA chr7:117,363,307-117,427,507(-) |
| chr7:3 | A00839:75:H7MFFDSXY:3:1378:6497:19006 | chr7:158,326,625-158,326,922 | nrna | 35M2I27M2D86M | 298 | 1 | 5 | 0.18277612328529358 | ambig_same_strand | unspliced | nRNA chr7:157,539,055-158,587,823(-) |
| chr7:3 | A00839:75:H7MFFDSXY:3:1666:13919:21277 | chr7:90,163,010-90,163,350 | nrna | 150M | -341 | 1 | 0 | 0.8226590156555176 | ambig_opp_strand | unspliced | nRNA chr7:90,161,265-90,164,829(+) |
| chr7:3 | A00839:75:H7MFFDSXY:3:2458:2130:31140 | chr7:128,493,573-128,493,874 | nrna | 150M | -302 | 1 | 0 | 0.7546987533569336 | ambig_same_strand | unspliced | nRNA chr7:128,476,728-128,502,126(+) |
| chr3:3 | A00839:75:H7MFFDSXY:3:2437:16893:33786 | chr3:47,343,515-47,343,751 | nrna | 150M | -237 | 1 | 0 | 0.8950293660163879 | ambig_same_strand | unspliced | nRNA chr3:47,282,943-47,344,577(+) |
| chr3:3 | A00839:75:H7MFFDSXY:3:2619:32958:33332 | chr3:53,775,974-53,776,112 | nrna | 139M | -139 | 1 | 1 | 0.3855988085269928 | ambig_same_strand | unspliced | nRNA chr3:53,743,011-53,792,062(+) |
| chr3:3 | A00839:75:H7MFFDSXY:3:1233:16206:23296 | chr3:53,800,253-53,800,507 | nrna | 150M | -255 | 1 | 0 | 0.9603784084320068 | ambig_same_strand | unspliced | nRNA chr3:53,799,890-53,813,712(+) |
| chr3:3 | A00839:75:H7MFFDSXY:3:2552:23764:25801 | chr3:47,891,475-47,891,772 | mrna | 150M | 298 | 1 | 1 | 0.8829472661018372 | ambig_same_strand | unspliced | mRNA MAP4 ENST00000497735.1 |
| chr4:3 | A00839:75:H7MFFDSXY:3:1168:30309:18067 | chr4:99,949,635-99,949,866 | nrna | 150M | 232 | 1 | 0 | 0.9068192839622498 | ambig_same_strand | unspliced | nRNA chr4:99,948,777-99,950,268(-) |
| chr4:3 | A00839:75:H7MFFDSXY:3:2346:21956:29465 | chr4:87,612,481-87,612,737 | mrna | 150M | -257 | 1 | 0 | 0.8435884714126587 | ambig_opp_strand | unspliced | mRNA DSPP ENST00000651931.1 |
| chr4:3 | A00839:75:H7MFFDSXY:3:1125:23475:34225 | chr4:67,825,627-67,825,993 | nrna | 150M | 367 | 1 | 0 | 0.4754561185836792 | ambig_opp_strand | unspliced | nRNA chr4:67,820,875-67,833,627(-) |
| chr4:3 | A00839:75:H7MFFDSXY:3:1472:7211:12414 | chr4:78,537,090-78,537,397 | nrna | 150M | -308 | 1 | 1 | 0.7230608463287354 | ambig_same_strand | unspliced | nRNA chr4:78,057,322-78,544,269(+) |
| chr5:3 | A00839:75:H7MFFDSXY:3:1134:16306:3114 | chr5:141,200,762-141,201,109 | nrna | 150M | 348 | 1 | 0 | 0.6548960208892822 | ambig_opp_strand | unspliced | nRNA chr5:141,136,682-141,241,744(-) |
| chr5:3 | A00839:75:H7MFFDSXY:3:1163:14262:25441 | chr5:141,183,091-141,183,317 | nrna | 150M | 227 | 1 | 0 | 0.9741189479827881 | ambig_opp_strand | unspliced | nRNA chr5:141,136,682-141,241,744(-) |
| chr5:3 | A00839:75:H7MFFDSXY:3:1512:7401:2409 | chr5:141,235,450-141,235,769 | nrna | 150M | -320 | 1 | 0 | 0.674885094165802 | ambig_opp_strand | unspliced | nRNA ENST00000524813.2 |
| chr5:3 | A00839:75:H7MFFDSXY:3:2438:24144:12837 | chr5:140,652,894-140,653,167 | nrna | 150M | -274 | 1 | 0 | 0.9476761221885681 | ambig_same_strand | unspliced | nRNA chr5:140,647,828-140,662,480(+) |
| chr9:3 | A00839:75:H7MFFDSXY:3:1143:16740:3208 | chr9:120,471,757-120,471,980 | nrna | 150M | 224 | 1 | 0 | 0.7599226236343384 | ambig_same_strand | unspliced | nRNA chr9:120,388,874-120,475,733(-) |
| chr9:3 | A00839:75:H7MFFDSXY:3:1647:25554:35289 | chr9:32,633,742-32,634,092 | nrna | 150M | 351 | 1 | 0 | 0.9758726954460144 | ambig_opp_strand | unspliced | nRNA ENST00000242310.4 |
| chr9:3 | A00839:75:H7MFFDSXY:3:2358:21423:33833 | chr9:110,481,103-110,481,314 | nrna | 150M | 212 | 1 | 0 | 0.768936812877655 | ambig_same_strand | unspliced | nRNA chr9:110,365,247-110,579,741(-) |
| chr9:3 | A00839:75:H7MFFDSXY:3:1368:17580:27993 | chr9:112,380,123-112,380,342 | nrna | 150M | -220 | 1 | 0 | 0.7718521356582642 | ambig_same_strand | unspliced | nRNA chr9:112,380,079-112,472,405(+) |
| chr8:3 | A00839:75:H7MFFDSXY:3:1517:4110:30279 | chr8:99,891,808-99,892,048 | nrna | 150M | 241 | 1 | 1 | 0.508952260017395 | ambig_same_strand | unspliced | nRNA chr8:99,874,995-99,893,707(-) |
| chr8:3 | A00839:75:H7MFFDSXY:3:1343:12310:18959 | chr8:72,937,822-72,937,992 | mrna | 150M | -171 | 1 | 0 | 0.9929520487785339 | ambig_same_strand | unspliced | mRNA KCNB2 ENST00000523207.2 |

## High-Confidence False-Positive Examples

These examples have the largest `ZW` winner weights among gDNA-source fragments called RNA.

| QNAME | Region | Detail | CIGAR | TLEN | NH | NM | ZW | ZC | ZS | Target |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A00839:75:H7MFFDSXY:3:2140:1479:20556 | chrX:31,126,583-79,061,231 | mrna | 150M | 47934649 | 1 | 2 | 1.0 | ambig_same_strand | unspliced | mRNA DMD ENST00000357033.9 |
| A00839:75:H7MFFDSXY:3:2164:12418:20087 | chr17:28,634,058-28,634,970 | mrna | 150M | -913 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA ENSG00000264044 ENST00000583787.1 |
| A00839:75:H7MFFDSXY:3:2557:32606:23985 | chr19:3,820,642-3,821,102 | nrna | 150M | 461 | 1 | 2 | 1.0 | ambig_same_strand | spliced_unannot | nRNA chr19:3,804,023-3,869,038(-) |
| A00839:75:H7MFFDSXY:3:1147:20772:10285 | chr1:165,708,138-165,708,350 | nrna | 7S85M58S | 213 | 2 | 7 | 1.0 | multimapper | spliced_unannot | nRNA chr1:165,698,749-165,709,968(+) |
| A00839:75:H7MFFDSXY:3:1219:25021:17613 | chr3:183,231,087-183,276,874 | mrna | 126M45659N3M | -45788 | 1 | 0 | 1.0 | ambig_same_strand | unspliced | mRNA MCF2L2 ENST00000488149.5 |
| A00839:75:H7MFFDSXY:3:2375:32425:32127 | chr17:21,689,501-21,689,717 | nrna | 150M | -217 | 2 | 0 | 1.0 | multimapper | spliced_unannot | nRNA chr17:21,376,356-21,419,870(+) |
| A00839:75:H7MFFDSXY:3:2611:18656:2988 | chr14:52,511,185-52,511,507 | mrna | 150M | 323 | 6 | 1 | 1.0 | multimapper | unspliced | mRNA TXNDC16 ENST00000281741.9 |
| A00839:75:H7MFFDSXY:3:2239:19253:36370 | chr17:41,084,030-41,084,372 | nrna | 150M | -343 | 2 | 2 | 1.0 | multimapper | spliced_unannot | nRNA ENST00000391417.6 |
| A00839:75:H7MFFDSXY:3:1628:29116:23328 | chr17:81,415,970-81,416,888 | nrna | 150M | -919 | 1 | 0 | 1.0 | ambig_same_strand | spliced_unannot | nRNA chr17:81,395,456-81,466,331(+) |
| A00839:75:H7MFFDSXY:3:2543:18846:25801 | chr20:35,000,864-35,001,318 | nrna | 35M22N115M | 455 | 1 | 0 | 1.0 | ambig_same_strand | spliced_unannot | nRNA chr20:34,955,867-35,002,437(+) |
| A00839:75:H7MFFDSXY:3:1324:23258:27931 | chr12:108,623,536-108,623,944 | mrna | 49M | 409 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA SELPLG ENST00000388962.4 |
| A00839:75:H7MFFDSXY:3:2674:17824:23938 | chr8:144,085,136-144,085,810 | mrna | 150M | -675 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA GPAA1 ENST00000703634.1 |
| A00839:75:H7MFFDSXY:3:2370:17562:36166 | chr12:8,917,751-8,921,612 | mrna | 41M1964N109M | 3862 | 1 | 1 | 1.0 | ambig_same_strand | unspliced | mRNA PHC1 ENST00000541181.1 |
| A00839:75:H7MFFDSXY:3:1268:6976:35587 | chr17:78,543,619-78,547,239 | nrna | 150M | 3621 | 2 | 3 | 1.0 | multimapper | spliced_unannot | nRNA chr17:78,423,696-78,577,396(-) |
| A00839:75:H7MFFDSXY:3:2155:22697:1219 | chr11:134,736,903-134,737,270 | nrna | 150M | 368 | 8 | 0 | 1.0 | multimapper | spliced_unannot | nRNA chr11:134,735,595-134,763,810(+) |
| A00839:75:H7MFFDSXY:3:2560:21314:22153 | chr19:54,219,600-54,239,734 | nrna | 76M19744N74M | -20135 | 1 | 0 | 1.0 | ambig_same_strand | spliced_unannot | nRNA chr19:54,200,857-54,249,003(+) |
| A00839:75:H7MFFDSXY:3:2174:6922:1892 | chr6:36,495,798-36,496,037 | mrna | 150M | 240 | 3 | 0 | 1.0 | multimapper | unspliced | mRNA STK38 ENST00000229812.8 |
| A00839:75:H7MFFDSXY:3:1163:25012:32816 | chr5:81,344,758-81,345,174 | nrna | 91M21N59M | -417 | 1 | 1 | 1.0 | ambig_same_strand | spliced_unannot | nRNA chr5:81,329,995-81,394,134(-) |
| A00839:75:H7MFFDSXY:3:1633:29315:12117 | chr1:31,756,668-31,759,323 | mrna | 148M2505N2M | 2656 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA ADGRB2 ENST00000398556.7 |
| A00839:75:H7MFFDSXY:3:1671:13078:27586 | chr8:58,425,881-58,426,519 | mrna | 150M | 639 | 1 | 1 | 1.0 | unambig | spliced_annot | mRNA PTPN11P2 ENST00000510501.1 |
| A00839:75:H7MFFDSXY:3:2226:30789:24189 | chr4:635,285-635,876 | nrna | 150M | -592 | 1 | 0 | 1.0 | ambig_same_strand | spliced_unannot | nRNA chr4:625,572-670,782(+) |
| A00839:75:H7MFFDSXY:3:2248:16721:9377 | chr6:32,178,601-32,180,054 | mrna | 49M113N55M77N46M | -1454 | 1 | 3 | 1.0 | unambig | spliced_annot | mRNA RNF5 ENST00000375094.4 |
| A00839:75:H7MFFDSXY:3:2265:15031:3787 | chr21:30,741,884-30,742,616 | nrna | 150M | -733 | 1 | 2 | 1.0 | unambig | spliced_unannot | nRNA ENST00000454921.1 |
| A00839:75:H7MFFDSXY:3:1138:22101:31344 | chr20:62,417,123-62,417,528 | nrna | 150M | 406 | 8 | 1 | 1.0 | multimapper | spliced_unannot | nRNA chr20:62,410,236-62,427,539(-) |
| A00839:75:H7MFFDSXY:3:1349:12255:3834 | chr3:48,468,215-48,469,656 | mrna | 150M | -1442 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA SHISA5 ENST00000466424.5 |
| A00839:75:H7MFFDSXY:3:1114:2682:17065 | chr21:46,334,350-46,334,707 | mrna | 52M3D1M75N97M | -358 | 1 | 4 | 1.0 | unambig | spliced_annot | mRNA PCNT ENST00000695526.1 |
| A00839:75:H7MFFDSXY:3:1633:25400:19366 | chr12:103,985,687-103,987,320 | mrna | 42M1219N108M | 1634 | 1 | 0 | 1.0 | ambig_same_strand | spliced_annot | mRNA TDG ENST00000392872.8 |
| A00839:75:H7MFFDSXY:3:1609:2230:29653 | chr20:62,338,610-62,338,794 | nrna | 150M | -185 | 2 | 1 | 1.0 | multimapper | spliced_unannot | nRNA chr20:62,309,058-62,367,312(-) |
| A00839:75:H7MFFDSXY:3:2531:23900:30514 | chr12:108,623,720-108,624,005 | mrna | 62M30N88M | -286 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA SELPLG ENST00000388962.4 |
| A00839:75:H7MFFDSXY:3:1113:2157:22670 | chr1:7,829,878-7,830,233 | mrna | 74M54N76M | 356 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA PER3 ENST00000614998.4 |

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

- Error categories: `results/vcap_gdna_false_rna_hotspots_2026-05-16/gdna_false_rna_categories.tsv`
- Hotspot windows: `results/vcap_gdna_false_rna_hotspots_2026-05-16/gdna_false_rna_windows.tsv`
- Hotspot EM loci: `results/vcap_gdna_false_rna_hotspots_2026-05-16/gdna_false_rna_loci.tsv`
- Assigned RNA targets: `results/vcap_gdna_false_rna_hotspots_2026-05-16/gdna_false_rna_targets.tsv`
- Representative fragments: `results/vcap_gdna_false_rna_hotspots_2026-05-16/sample_false_rna_fragments.tsv`
- High-confidence examples: `results/vcap_gdna_false_rna_hotspots_2026-05-16/high_confidence_false_rna_fragments.tsv`
- Summary JSON: `results/vcap_gdna_false_rna_hotspots_2026-05-16/summary.json`
