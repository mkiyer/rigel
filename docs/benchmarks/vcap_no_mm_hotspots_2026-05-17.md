# VCaP gDNA False-RNA Hotspot Analysis - 2026-05-16

BAM: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/no_mm/annotated.bam`

Truth source is derived from query-name flowcell ID. This analysis focuses on gDNA-source fragments (`H7MFFDSXY`) that Rigel assigned to an RNA pool (`mRNA` or `nRNA`). Counting rule is one primary read1 record per fragment.

Fragments counted: 32,218,601

gDNA-source fragments: 18,228,677

gDNA -> RNA false positives: 2,640,324 (14.48% of gDNA-source fragments)

False-positive split: 1,310,889 mRNA (49.65% of false RNA) and 1,329,435 nRNA (50.35% of false RNA).

## What Dominates

The false-positive RNA problem is concentrated in unspliced, genic ambiguous fragments. Across false RNA calls, 97.10% are `ZS=unspliced` and 98.92% are either `ambig_same_strand` or `ambig_opp_strand`.

The top 10 EM loci account for 3.58% of all gDNA -> RNA false positives. The top 10 genomic 10,000 bp windows account for 1.09%. So this is not only a few pathological coordinates; it is a recurrent ambiguous-region failure mode, with some strong hotspots.

## Error Categories

| Pred detail | ZC | ZS | Count | False RNA frac | gDNA source frac |
| --- | --- | --- | --- | --- | --- |
| mrna | ambig_same_strand | unspliced | 1,043,881 | 39.54% | 5.73% |
| nrna | ambig_same_strand | unspliced | 862,648 | 32.67% | 4.73% |
| nrna | ambig_opp_strand | unspliced | 388,357 | 14.71% | 2.13% |
| mrna | ambig_opp_strand | unspliced | 244,051 | 9.24% | 1.34% |
| nrna | ambig_same_strand | spliced_implicit | 41,915 | 1.59% | 0.23% |
| nrna | unambig | unspliced | 24,897 | 0.94% | 0.14% |
| mrna | ambig_same_strand | spliced_annot | 9,159 | 0.35% | 0.05% |
| nrna | ambig_opp_strand | spliced_implicit | 6,793 | 0.26% | 0.04% |
| mrna | ambig_same_strand | spliced_implicit | 6,604 | 0.25% | 0.04% |
| nrna | ambig_same_strand | spliced_unannot | 4,137 | 0.16% | 0.02% |
| mrna | unambig | spliced_annot | 3,506 | 0.13% | 0.02% |
| mrna | ambig_same_strand | spliced_unannot | 1,820 | 0.07% | 0.01% |
| mrna | ambig_opp_strand | spliced_implicit | 1,697 | 0.06% | 0.01% |
| nrna | ambig_opp_strand | spliced_unannot | 535 | 0.02% | 0.00% |
| mrna | ambig_opp_strand | spliced_unannot | 171 | 0.01% | 0.00% |
| nrna | unambig | spliced_unannot | 153 | 0.01% | 0.00% |

## Top 10kb Genomic Windows

| Region | False RNA | False rate | gDNA total | mRNA | nRNA | Top target | Top category |
| --- | --- | --- | --- | --- | --- | --- | --- |
| chr11:65,500,001-65,510,000 | 3,822 | 90.87% | 4,206 | 1,981 | 1,841 | nRNA ENST00000698129.1 | mrna/ambig_opp_strand/unspliced |
| chrX:67,540,001-67,550,000 | 3,270 | 49.66% | 6,585 | 3,256 | 14 | mRNA AR ENST00000396043.4 | mrna/ambig_same_strand/unspliced |
| chr11:64,220,001-64,230,000 | 2,994 | 51.54% | 5,809 | 762 | 2,232 | nRNA chr11:64,218,896-64,223,857(+) | nrna/ambig_same_strand/unspliced |
| chr1:152,350,001-152,360,000 | 2,868 | 62.21% | 4,610 | 1,503 | 1,365 | mRNA FLG2 ENST00000388718.5 | mrna/ambig_opp_strand/unspliced |
| chr9:106,920,001-106,930,000 | 2,839 | 59.17% | 4,798 | 2,689 | 150 | mRNA ZNF462 ENST00000277225.10 | mrna/ambig_same_strand/unspliced |
| chr2:178,560,001-178,570,000 | 2,825 | 71.54% | 3,949 | 249 | 2,576 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr2:178,570,001-178,580,000 | 2,641 | 71.05% | 3,717 | 197 | 2,444 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr2:178,530,001-178,540,000 | 2,627 | 71.62% | 3,668 | 351 | 2,276 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr11:64,230,001-64,240,000 | 2,542 | 55.06% | 4,617 | 186 | 2,356 | nRNA chr11:64,229,213-64,234,292(-) | nrna/ambig_opp_strand/unspliced |
| chr2:178,540,001-178,550,000 | 2,291 | 72.78% | 3,148 | 116 | 2,175 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr21:33,550,001-33,560,000 | 1,986 | 48.18% | 4,122 | 1,885 | 101 | mRNA SON ENST00000300278.8 | mrna/ambig_same_strand/unspliced |
| chr11:64,110,001-64,120,000 | 1,903 | 46.85% | 4,062 | 1,408 | 495 | mRNA FLRT1 ENST00000246841.3 | mrna/ambig_opp_strand/unspliced |
| chr11:66,330,001-66,340,000 | 1,841 | 57.41% | 3,207 | 77 | 1,764 | nRNA chr11:66,334,493-66,339,875(+) | nrna/ambig_opp_strand/unspliced |
| chr11:64,260,001-64,270,000 | 1,834 | 34.04% | 5,387 | 66 | 1,768 | nRNA chr11:64,251,529-64,267,923(+) | nrna/ambig_same_strand/unspliced |
| chr2:178,610,001-178,620,000 | 1,827 | 68.79% | 2,656 | 15 | 1,812 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr9:76,700,001-76,710,000 | 1,815 | 81.79% | 2,219 | 1,006 | 809 | nRNA chr9:76,704,847-76,710,661(+) | mrna/ambig_opp_strand/unspliced |
| chr11:66,850,001-66,860,000 | 1,729 | 48.35% | 3,576 | 669 | 1,060 | nRNA chr11:66,857,882-66,860,475(-) | nrna/ambig_opp_strand/unspliced |
| chr2:178,590,001-178,600,000 | 1,716 | 68.20% | 2,516 | 56 | 1,660 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr2:178,550,001-178,560,000 | 1,692 | 70.62% | 2,396 | 77 | 1,615 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr13:109,780,001-109,790,000 | 1,653 | 42.70% | 3,871 | 1,653 | 0 | mRNA IRS2 ENST00000375856.5 | mrna/ambig_same_strand/unspliced |
| chr2:178,600,001-178,610,000 | 1,536 | 67.22% | 2,285 | 25 | 1,511 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr11:64,250,001-64,260,000 | 1,513 | 33.21% | 4,556 | 17 | 1,496 | nRNA chr11:64,251,529-64,267,923(+) | nrna/ambig_same_strand/unspliced |
| chr2:178,580,001-178,590,000 | 1,509 | 65.92% | 2,289 | 49 | 1,460 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr17:29,570,001-29,580,000 | 1,504 | 51.33% | 2,930 | 153 | 1,351 | nRNA chr17:29,560,546-29,707,090(+) | nrna/ambig_opp_strand/unspliced |
| chr2:178,770,001-178,780,000 | 1,486 | 66.91% | 2,221 | 47 | 1,439 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr8:76,850,001-76,860,000 | 1,482 | 50.75% | 2,920 | 1,457 | 25 | mRNA ZFHX4 ENST00000518282.5 | mrna/ambig_same_strand/unspliced |
| chr11:63,760,001-63,770,000 | 1,469 | 45.85% | 3,204 | 121 | 1,348 | nRNA chr11:63,759,891-63,768,775(-) | nrna/ambig_same_strand/unspliced |
| chr11:67,250,001-67,260,000 | 1,467 | 55.57% | 2,640 | 501 | 966 | nRNA chr11:67,252,335-67,261,545(-) | nrna/ambig_opp_strand/unspliced |
| chr13:110,500,001-110,510,000 | 1,428 | 56.67% | 2,520 | 17 | 1,411 | nRNA chr13:110,502,574-110,508,179(-) | nrna/ambig_opp_strand/unspliced |
| chr2:185,790,001-185,800,000 | 1,410 | 67.11% | 2,101 | 716 | 694 | mRNA FSIP2 ENST00000424728.6 | mrna/ambig_opp_strand/unspliced |

## Local RNA Evidence Stratification

False RNA calls are not confined to windows with abundant RNA-source reads. This matters
because it separates two mechanisms: real local RNA expression can pull ambiguous gDNA
into RNA, but ambiguous gDNA can also self-seed RNA components in windows with little or
no RNA evidence.

| RNA-source fragments in 10kb window | Windows | False RNA | False RNA frac | Weighted false rate |
| --- | --- | --- | --- | --- |
| 0 | 15,468 | 333,869 | 12.65% | 13.40% |
| 1-100 | 25,737 | 777,442 | 29.44% | 14.17% |
| 101-1,000 | 16,222 | 975,572 | 36.95% | 17.38% |
| >1,000 | 2,958 | 553,441 | 20.96% | 23.69% |

## Top EM Loci

The region span is the observed fragment span among false-positive fragments in the locus, not the full locus extent. Components such as `chr2:3`, `chr1:3`, and `chr7:3` are chromosome-scale mega-components rather than precise biological loci, so the 10kb windows and target tables above are the sharper coordinates for assay triage.

| Locus | Observed region | False RNA | False frac | mRNA | nRNA | Same | Opp | Mean ZW | ZW >= .90 | Top target |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr2:12994 | chr2:24,586,848-186,333,274 | 32,893 | 1.25% | 1,704 | 31,189 | 668 | 32,224 | 0.519 | 13.46% | nRNA chr2:178,525,988-178,804,642(-) |
| chr11:3419 | chr11:644,233-114,450,251 | 14,507 | 0.55% | 4,801 | 9,706 | 5,609 | 8,871 | 0.567 | 22.14% | nRNA ENST00000698129.1 |
| chr13:6465 | chr13:96,090,628-110,512,596 | 9,744 | 0.37% | 450 | 9,294 | 6,778 | 2,965 | 0.554 | 15.95% | nRNA chr13:110,149,060-110,299,118(-) |
| chr2:12263 | chr2:3,651,440-240,805,352 | 6,514 | 0.25% | 2,112 | 4,402 | 5,032 | 1,428 | 0.503 | 10.41% | nRNA chr2:169,127,108-169,362,534(-) |
| chr17:9590 | chr17:29,562,177-30,185,966 | 6,196 | 0.23% | 1,685 | 4,511 | 1,595 | 4,600 | 0.590 | 13.51% | nRNA chr17:29,560,546-29,707,090(+) |
| chr5:18408 | chr5:107,380,994-141,248,274 | 5,981 | 0.23% | 135 | 5,846 | 216 | 5,659 | 0.752 | 37.99% | nRNA chr5:141,136,682-141,241,744(-) |
| chr1:1437 | chr1:152,124,320-152,412,433 | 5,574 | 0.21% | 3,677 | 1,897 | 836 | 4,738 | 0.678 | 9.35% | mRNA FLG2 ENST00000388718.5 |
| chr11:4278 | chr11:67,351,592-85,758,001 | 4,584 | 0.17% | 1,175 | 3,409 | 2,481 | 2,096 | 0.505 | 2.42% | nRNA chr11:67,400,854-67,421,183(-) |
| chr5:18628 | chr5:108,713,551-141,512,199 | 4,515 | 0.17% | 372 | 4,143 | 4,367 | 148 | 0.646 | 3.06% | nRNA ENST00000611950.1 |
| chr17:9377 | chr17:10,300,785-46,371,505 | 4,128 | 0.16% | 141 | 3,987 | 19 | 4,105 | 0.526 | 4.60% | nRNA chr17:10,383,155-10,537,862(+) |
| chrX:23984 | chrX:67,521,831-67,853,372 | 4,117 | 0.16% | 3,964 | 153 | 4,116 | 1 | 0.365 | 9.79% | mRNA AR ENST00000396043.4 |
| chr9:22928 | chr9:21,069,650-123,930,130 | 4,107 | 0.16% | 2,255 | 1,852 | 1,890 | 2,217 | 0.642 | 29.90% | nRNA chr9:76,704,847-76,710,661(+) |
| chr17:9572 | chr17:28,473,423-28,645,234 | 4,032 | 0.15% | 895 | 3,137 | 2,169 | 1,861 | 0.535 | 5.78% | nRNA chr17:28,614,881-28,645,454(-) |
| chr9:22644 | chr9:27,723,901-107,012,899 | 3,855 | 0.15% | 3,312 | 543 | 3,657 | 198 | 0.614 | 31.60% | mRNA ZNF462 ENST00000277225.10 |
| chr21:14756 | chr21:11,628,714-33,915,950 | 3,853 | 0.15% | 2,504 | 1,349 | 3,010 | 843 | 0.466 | 5.76% | mRNA SON ENST00000300278.8 |
| chr3:15589 | chr3:9,749,983-195,443,010 | 3,828 | 0.14% | 1,766 | 2,062 | 1,820 | 2,004 | 0.548 | 15.44% | mRNA IGSF10 ENST00000282466.4 |
| chr11:3535 | chr11:5,248,270-62,984,861 | 3,735 | 0.14% | 1,293 | 2,442 | 1,499 | 2,235 | 0.722 | 38.47% | nRNA ENST00000380184.2 |
| chr12:5463 | chr12:57,096,217-84,580,459 | 3,623 | 0.14% | 320 | 3,303 | 3,355 | 266 | 0.609 | 3.81% | nRNA chr12:57,128,482-57,213,361(+) |
| chr11:4167 | chr11:63,974,619-64,166,104 | 3,517 | 0.13% | 1,756 | 1,761 | 1,157 | 2,359 | 0.688 | 31.11% | mRNA FLRT1 ENST00000246841.3 |
| chr11:4067 | chr11:59,712,942-64,269,089 | 3,433 | 0.13% | 85 | 3,348 | 3,345 | 47 | 0.801 | 43.98% | nRNA chr11:64,251,529-64,267,923(+) |
| chr5:18625 | chr5:140,786,133-141,011,127 | 3,126 | 0.12% | 659 | 2,467 | 1,697 | 1,429 | 0.644 | 9.72% | nRNA ENST00000409700.4 |
| chr6:19130 | chr6:26,026,894-166,439,391 | 3,064 | 0.12% | 1,038 | 2,026 | 2,193 | 829 | 0.609 | 24.87% | nRNA chr6:26,457,949-26,469,653(+) |
| chr11:3698 | chr11:26,249,163-94,492,977 | 2,990 | 0.11% | 867 | 2,123 | 2,793 | 196 | 0.563 | 1.71% | nRNA chr11:68,047,283-68,050,880(+) |
| chr1:254 | chr1:13,893,297-228,495,719 | 2,891 | 0.11% | 1,419 | 1,472 | 2,751 | 83 | 0.527 | 12.80% | nRNA chr1:154,582,056-154,599,446(-) |
| chr3:15741 | chr3:38,039,218-177,106,114 | 2,800 | 0.11% | 507 | 2,293 | 1,476 | 1,319 | 0.574 | 10.82% | nRNA chr3:38,039,204-38,124,025(+) |
| chr2:13468 | chr2:185,738,484-185,833,291 | 2,776 | 0.11% | 1,919 | 857 | 632 | 2,144 | 0.813 | 47.73% | mRNA FSIP2 ENST00000424728.6 |
| chr9:23236 | chr9:114,862,696-115,401,293 | 2,772 | 0.10% | 754 | 2,018 | 20 | 2,752 | 0.499 | 0.87% | nRNA chr9:114,850,967-115,121,030(+) |
| chr3:15939 | chr3:52,495,346-190,656,792 | 2,755 | 0.10% | 328 | 2,427 | 2,616 | 135 | 0.509 | 4.21% | nRNA chr3:52,495,337-52,524,495(+) |
| chr7:20901 | chr7:82,821,985-143,352,083 | 2,714 | 0.10% | 2,244 | 470 | 2,712 | 1 | 0.677 | 40.38% | mRNA PCLO ENST00000423517.6 |
| chr8:22083 | chr8:76,603,277-92,818,668 | 2,688 | 0.10% | 2,608 | 80 | 2,686 | 2 | 0.694 | 44.31% | mRNA ZFHX4 ENST00000518282.5 |

## Top Assigned RNA Targets

| Target | Observed region | False RNA | False frac | mRNA | nRNA | Same | Opp | Mean ZW |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| nRNA chr2:178,525,988-178,804,642(-) | chr2:178,526,839-178,804,618 | 14,782 | 0.56% | 0 | 14,782 | 440 | 14,342 | 0.757 |
| nRNA chr2:178,523,417-178,776,534(+) | chr2:178,526,721-178,776,531 | 7,662 | 0.29% | 0 | 7,662 | 0 | 7,662 | 0.507 |
| nRNA chr11:64,251,529-64,267,923(+) | chr11:64,251,539-64,267,707 | 3,264 | 0.12% | 0 | 3,264 | 3,264 | 0 | 0.809 |
| nRNA chr13:110,149,060-110,299,118(-) | chr13:110,149,064-110,299,077 | 2,492 | 0.09% | 0 | 2,492 | 2,492 | 0 | 0.597 |
| mRNA ZFHX4 ENST00000518282.5 | chr8:76,704,057-92,818,668 | 2,184 | 0.08% | 2,184 | 0 | 2,183 | 1 | 0.787 |
| mRNA ZNF462 ENST00000277225.10 | chr9:29,052,401-107,012,899 | 2,176 | 0.08% | 2,176 | 0 | 2,151 | 25 | 0.734 |
| mRNA AHNAK ENST00000378024.9 | chr11:62,516,540-62,534,074 | 2,130 | 0.08% | 2,130 | 0 | 2,130 | 0 | 0.842 |
| nRNA chr12:57,128,482-57,213,361(+) | chr12:57,128,772-57,212,660 | 2,094 | 0.08% | 0 | 2,094 | 2,009 | 85 | 0.632 |
| nRNA chr17:29,560,546-29,707,090(+) | chr17:29,562,309-29,703,260 | 2,088 | 0.08% | 0 | 2,088 | 250 | 1,838 | 0.575 |
| nRNA chr2:178,528,764-178,620,148(+) | chr2:178,528,774-178,619,987 | 2,065 | 0.08% | 0 | 2,065 | 0 | 2,065 | 0.180 |
| nRNA chr11:64,823,051-64,844,653(-) | chr11:64,687,206-64,844,644 | 1,968 | 0.07% | 0 | 1,968 | 1,968 | 0 | 0.795 |
| mRNA MUC16 ENST00000397910.8 | chr19:8,934,900-8,981,296 | 1,888 | 0.07% | 1,888 | 0 | 1,888 | 0 | 0.554 |
| mRNA PCLO ENST00000423517.6 | chr7:82,821,985-83,162,880 | 1,888 | 0.07% | 1,888 | 0 | 1,888 | 0 | 0.787 |
| nRNA ENST00000698129.1 | chr11:65,499,313-65,506,876 | 1,824 | 0.07% | 0 | 1,824 | 0 | 1,824 | 0.864 |
| nRNA chr5:141,136,682-141,241,744(-) | chr5:141,136,747-141,241,485 | 1,804 | 0.07% | 0 | 1,804 | 3 | 1,801 | 0.678 |
| mRNA AR ENST00000396043.4 | chrX:67,521,831-67,724,284 | 1,769 | 0.07% | 1,769 | 0 | 1,769 | 0 | 0.479 |
| mRNA FSIP2 ENST00000424728.6 | chr2:185,788,659-185,833,291 | 1,730 | 0.07% | 1,730 | 0 | 631 | 1,099 | 0.822 |
| mRNA BSN ENST00000296452.5 | chr3:49,624,975-49,671,352 | 1,723 | 0.07% | 1,723 | 0 | 1,722 | 0 | 0.847 |
| mRNA IRS2 ENST00000375856.5 | chr13:109,755,789-109,786,266 | 1,705 | 0.06% | 1,705 | 0 | 1,491 | 214 | 0.877 |
| mRNA FAT3 ENST00000409404.6 | chr11:92,352,096-92,891,411 | 1,697 | 0.06% | 1,697 | 0 | 1,683 | 14 | 0.841 |
| mRNA FAT1 ENST00000441802.7 | chr4:83,052,074-186,709,837 | 1,572 | 0.06% | 1,572 | 0 | 1,568 | 1 | 0.847 |
| mRNA FLG2 ENST00000388718.5 | chr1:152,350,343-152,357,635 | 1,503 | 0.06% | 1,503 | 0 | 0 | 1,503 | 0.767 |
| mRNA SON ENST00000300278.8 | chr21:33,546,227-33,560,571 | 1,467 | 0.06% | 1,467 | 0 | 1,404 | 63 | 0.606 |
| mRNA FAT4 ENST00000394329.9 | chr4:125,315,073-125,492,011 | 1,434 | 0.05% | 1,434 | 0 | 1,434 | 0 | 0.798 |
| mRNA FLRT1 ENST00000246841.3 | chr11:64,103,236-64,118,914 | 1,409 | 0.05% | 1,409 | 0 | 0 | 1,409 | 0.851 |
| nRNA chr13:110,477,922-110,513,027(+) | chr13:110,477,931-110,512,596 | 1,395 | 0.05% | 0 | 1,395 | 972 | 423 | 0.613 |
| nRNA chr11:64,340,203-64,357,534(+) | chr11:64,340,205-64,383,523 | 1,395 | 0.05% | 0 | 1,395 | 1,389 | 6 | 0.534 |
| nRNA chr11:64,229,213-64,234,292(-) | chr11:64,229,224-64,234,291 | 1,379 | 0.05% | 0 | 1,379 | 0 | 1,379 | 0.804 |
| mRNA SACS ENST00000402364.1 | chr13:23,249,827-23,358,337 | 1,303 | 0.05% | 1,303 | 0 | 1,302 | 1 | 0.642 |
| mRNA ZFP62 ENST00000502412.2 | chr5:107,995,110-180,851,489 | 1,250 | 0.05% | 1,250 | 0 | 1,249 | 1 | 0.823 |

## Representative False-Positive Fragments From Top Loci

| Locus | QNAME | Region | Detail | CIGAR | TLEN | NH | NM | ZW | ZC | ZS | Target |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr2:12994 | A00839:75:H7MFFDSXY:3:1552:32497:13463 | chr2:178,718,006-178,718,262 | nrna | 150M | -257 | 1 | 0 | 0.803736686706543 | ambig_opp_strand | unspliced | nRNA chr2:178,523,417-178,776,534(+) |
| chr2:12994 | A00839:75:H7MFFDSXY:3:1550:14687:8484 | chr2:178,534,544-178,534,787 | nrna | 150M | 244 | 1 | 0 | 0.742862343788147 | ambig_opp_strand | unspliced | nRNA chr2:178,525,988-178,804,642(-) |
| chr2:12994 | A00839:75:H7MFFDSXY:3:2659:6225:25175 | chr2:178,541,313-178,541,621 | nrna | 150M | -309 | 1 | 0 | 0.23120346665382385 | ambig_opp_strand | unspliced | nRNA chr2:178,523,417-178,776,534(+) |
| chr2:12994 | A00839:75:H7MFFDSXY:3:1237:23294:34507 | chr2:178,568,782-178,569,136 | nrna | 150M | -355 | 1 | 0 | 0.11680469661951065 | ambig_opp_strand | unspliced | nRNA chr2:178,540,115-178,582,663(+) |
| chr11:3419 | A00839:75:H7MFFDSXY:3:2224:27778:12054 | chr11:65,223,531-65,223,708 | nrna | 150M | -178 | 1 | 0 | 0.0813809484243393 | ambig_same_strand | unspliced | nRNA chr11:65,213,839-65,226,040(+) |
| chr11:3419 | A00839:75:H7MFFDSXY:3:1234:31620:12477 | chr11:64,219,947-64,220,282 | nrna | 150M | -336 | 1 | 0 | 0.39476507902145386 | ambig_same_strand | unspliced | nRNA chr11:64,219,509-64,220,669(+) |
| chr11:3419 | A00839:75:H7MFFDSXY:3:2675:25003:35055 | chr11:64,225,500-64,225,802 | nrna | 150M | 303 | 1 | 0 | 0.45932137966156006 | ambig_same_strand | unspliced | nRNA chr11:64,223,798-64,226,158(-) |
| chr11:3419 | A00839:75:H7MFFDSXY:3:1624:11686:30154 | chr11:65,500,844-65,501,189 | mrna | 150M | -346 | 1 | 0 | 0.7107172608375549 | ambig_opp_strand | unspliced | mRNA MALAT1 ENST00000710933.1 |
| chr13:6465 | A00839:75:H7MFFDSXY:3:1507:24306:15593 | chr13:110,198,547-110,198,700 | nrna | 150M | 154 | 1 | 0 | 0.6584106683731079 | ambig_same_strand | unspliced | nRNA chr13:110,149,060-110,299,118(-) |
| chr13:6465 | A00839:75:H7MFFDSXY:3:1414:8070:35759 | chr13:110,457,297-110,457,573 | nrna | 150M | -277 | 1 | 0 | 0.27382639050483704 | ambig_opp_strand | unspliced | nRNA chr13:110,307,283-110,504,595(+) |
| chr13:6465 | A00839:75:H7MFFDSXY:3:2525:22281:28933 | chr13:110,511,995-110,512,162 | nrna | 150M | -168 | 1 | 0 | 0.5444831252098083 | ambig_same_strand | unspliced | nRNA ENST00000648222.1 |
| chr13:6465 | A00839:75:H7MFFDSXY:3:2561:30933:34053 | chr13:110,163,363-110,163,584 | nrna | 150M | 222 | 1 | 0 | 0.6116864681243896 | ambig_same_strand | unspliced | nRNA chr13:110,149,060-110,299,118(-) |
| chr2:12263 | A00839:75:H7MFFDSXY:3:1516:25446:1188 | chr2:48,582,035-48,582,378 | mrna | 150M | -344 | 1 | 0 | 0.13352765142917633 | ambig_same_strand | unspliced | mRNA STON1 ENST00000404752.6 |
| chr2:12263 | A00839:75:H7MFFDSXY:3:2627:2944:30232 | chr2:201,271,391-201,271,601 | nrna | 150M | -211 | 1 | 0 | 0.07111163437366486 | ambig_same_strand | unspliced | nRNA chr2:201,261,951-201,272,776(+) |
| chr2:12263 | A00839:75:H7MFFDSXY:3:2266:14832:2910 | chr2:48,646,676-48,647,001 | mrna | 150M | -326 | 1 | 0 | 0.8795478940010071 | ambig_opp_strand | unspliced | mRNA STON1-GTF2A1L ENST00000394751.5 |
| chr2:12263 | A00839:75:H7MFFDSXY:3:1678:25527:24001 | chr2:127,286,839-127,287,146 | nrna | 150M | 308 | 1 | 0 | 0.6087235808372498 | ambig_same_strand | unspliced | nRNA chr2:127,286,748-127,294,107(-) |
| chr17:9590 | A00839:75:H7MFFDSXY:3:1105:2284:25363 | chr17:29,636,022-29,636,305 | mrna | 150M | 284 | 1 | 0 | 0.9675247669219971 | ambig_opp_strand | unspliced | mRNA SSH2 ENST00000649863.1 |
| chr17:9590 | A00839:75:H7MFFDSXY:3:1351:5972:23766 | chr17:29,612,073-29,612,410 | nrna | 150M | -338 | 1 | 0 | 0.23531405627727509 | ambig_same_strand | unspliced | nRNA chr17:29,560,546-29,707,090(+) |
| chr17:9590 | A00839:75:H7MFFDSXY:3:1506:25048:9518 | chr17:29,631,500-29,631,681 | nrna | 150M | -182 | 1 | 0 | 0.7121607065200806 | ambig_opp_strand | unspliced | nRNA chr17:29,560,546-29,707,090(+) |
| chr17:9590 | A00839:75:H7MFFDSXY:3:2645:30472:28620 | chr17:29,574,841-29,575,137 | nrna | 150M | -297 | 1 | 0 | 0.5534963011741638 | ambig_opp_strand | unspliced | nRNA chr17:29,560,546-29,707,090(+) |
| chr5:18408 | A00839:75:H7MFFDSXY:3:1163:14262:25441 | chr5:141,183,091-141,183,317 | nrna | 150M | 227 | 1 | 0 | 0.8757743835449219 | ambig_opp_strand | unspliced | nRNA chr5:141,136,682-141,241,744(-) |
| chr5:18408 | A00839:75:H7MFFDSXY:3:2258:19117:23171 | chr5:141,186,893-141,187,403 | nrna | 150M | 511 | 1 | 0 | 0.015096351504325867 | ambig_opp_strand | unspliced | nRNA chr5:141,183,400-141,201,396(-) |
| chr5:18408 | A00839:75:H7MFFDSXY:3:1512:7401:2409 | chr5:141,235,450-141,235,769 | nrna | 150M | -320 | 1 | 0 | 0.6418857574462891 | ambig_opp_strand | unspliced | nRNA ENST00000524813.2 |
| chr5:18408 | A00839:75:H7MFFDSXY:3:2462:6542:9909 | chr5:141,209,640-141,210,038 | nrna | 150M | -399 | 1 | 1 | 0.2515006363391876 | ambig_opp_strand | unspliced | nRNA chr5:141,208,696-141,211,537(+) |
| chr1:1437 | A00839:75:H7MFFDSXY:3:1436:25536:12712 | chr1:152,350,515-152,350,871 | nrna | 150M | 357 | 1 | 0 | 0.06421361863613129 | ambig_opp_strand | unspliced | nRNA chr1:152,348,734-152,360,006(-) |
| chr1:1437 | A00839:75:H7MFFDSXY:3:2130:1289:21762 | chr1:152,353,387-152,353,563 | nrna | 150M | -177 | 1 | 1 | 0.798764705657959 | ambig_opp_strand | unspliced | nRNA chr1:152,335,378-152,364,080(+) |
| chr1:1437 | A00839:75:H7MFFDSXY:3:2167:21920:30311 | chr1:152,155,615-152,155,895 | mrna | 150M | 281 | 1 | 1 | 0.8730872273445129 | ambig_same_strand | unspliced | mRNA RPTN ENST00000316073.3 |
| chr1:1437 | A00839:75:H7MFFDSXY:3:1208:28393:8735 | chr1:152,219,117-152,219,375 | mrna | 150M | 259 | 1 | 1 | 0.7690253853797913 | ambig_opp_strand | unspliced | mRNA HRNR ENST00000368801.4 |
| chr11:4278 | A00839:75:H7MFFDSXY:3:1352:22507:4178 | chr11:67,405,409-67,405,703 | nrna | 150M | -295 | 1 | 0 | 0.1007670909166336 | ambig_opp_strand | unspliced | nRNA chr11:67,405,247-67,410,089(+) |
| chr11:4278 | A00839:75:H7MFFDSXY:3:2111:13557:26224 | chr11:67,397,370-67,397,631 | nrna | 150M | -262 | 1 | 0 | 0.934227466583252 | ambig_same_strand | unspliced | nRNA chr11:67,395,669-67,398,410(+) |

## High-Confidence False-Positive Examples

These examples have the largest `ZW` winner weights among gDNA-source fragments called RNA.

| QNAME | Region | Detail | CIGAR | TLEN | NH | NM | ZW | ZC | ZS | Target |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A00839:75:H7MFFDSXY:3:2531:9372:8672 | chr19:4,538,574-4,539,090 | mrna | 150M | 517 | 1 | 1 | 1.0 | unambig | spliced_annot | mRNA LRG1 ENST00000586883.1 |
| A00839:75:H7MFFDSXY:3:2453:22019:36714 | chr1:7,829,834-7,830,156 | mrna | 118M54N32M | 323 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA PER3 ENST00000614998.4 |
| A00839:75:H7MFFDSXY:3:2276:23095:5321 | chr11:639,543-639,919 | mrna | 3M102N147M | 377 | 1 | 0 | 1.0 | ambig_same_strand | spliced_annot | mRNA DRD4 ENST00000176183.6 |
| A00839:75:H7MFFDSXY:3:1571:17056:14747 | chr22:50,530,318-50,530,945 | nrna | 150M | 628 | 1 | 1 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chr22:50,529,709-50,531,822(-) |
| A00839:75:H7MFFDSXY:3:2465:19316:4069 | chr22:22,635,852-22,636,197 | nrna | 150M | 346 | 1 | 9 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chr22:22,630,064-22,636,153(+) |
| A00839:75:H7MFFDSXY:3:1155:12057:5243 | chr5:169,882,852-169,883,168 | mrna | 126M24S | -317 | 1 | 1 | 1.0 | ambig_opp_strand | spliced_implicit | mRNA INSYN2B ENST00000377365.4 |
| A00839:75:H7MFFDSXY:3:1102:27697:34021 | chr14:24,238,329-24,238,961 | mrna | 150M | 633 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA GMPR2 ENST00000559479.1 |
| A00839:75:H7MFFDSXY:3:2249:26078:18223 | chr12:56,432,248-56,433,028 | mrna | 150M | 781 | 1 | 1 | 1.0 | unambig | spliced_annot | mRNA TIMELESS ENST00000229201.4 |
| A00839:75:H7MFFDSXY:3:1260:9498:7075 | chr17:81,700,422-81,700,713 | mrna | 110M178N4M | -292 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA HGS ENST00000677706.1 |
| A00839:75:H7MFFDSXY:3:1374:22833:3145 | chr20:35,000,854-35,001,279 | nrna | 150M | -426 | 1 | 0 | 1.0 | ambig_same_strand | spliced_unannot | nRNA chr20:34,955,867-35,002,437(+) |
| A00839:75:H7MFFDSXY:3:1648:14678:28855 | chr3:28,975,635-29,200,058 | mrna | 143M16922N7M | -224424 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA RBMS3 ENST00000636582.1 |
| A00839:75:H7MFFDSXY:3:1106:31195:22795 | chr8:144,057,560-144,058,049 | nrna | 150M | 490 | 1 | 0 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chr8:144,051,265-144,060,692(-) |
| A00839:75:H7MFFDSXY:3:1521:2754:27806 | chr3:93,470,611-113,420,085 | mrna | 150M | 19949475 | 1 | 9 | 1.0 | ambig_opp_strand | spliced_implicit | mRNA CFAP44 ENST00000393845.9 |
| A00839:75:H7MFFDSXY:3:1164:22607:5102 | chr6:33,417,861-33,418,272 | nrna | 150M | 412 | 1 | 0 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chr6:33,416,441-33,418,317(-) |
| A00839:75:H7MFFDSXY:3:2462:32199:29324 | chr6:70,149,718-70,150,217 | nrna | 150M | 500 | 1 | 0 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chr6:69,866,555-70,212,468(+) |
| A00839:75:H7MFFDSXY:3:1429:30255:16000 | chr3:196,782,753-196,806,610 | mrna | 86M2187N32M1195N32M | -23858 | 1 | 6 | 1.0 | unambig | spliced_annot | mRNA PAK2 ENST00000327134.7 |
| A00839:75:H7MFFDSXY:3:1477:29631:27759 | chr6:10,393,332-10,393,600 | mrna | 150M | 269 | 1 | 0 | 1.0 | ambig_same_strand | spliced_implicit | mRNA TFAP2A ENST00000461628.5 |
| A00839:75:H7MFFDSXY:3:1131:2193:24580 | chr17:36,599,949-36,600,308 | mrna | 150M | 360 | 1 | 0 | 1.0 | ambig_same_strand | spliced_annot | mRNA DHRS11 ENST00000617959.1 |
| A00839:75:H7MFFDSXY:3:1145:2320:1344 | chr17:47,819,951-47,820,370 | mrna | 150M | -420 | 1 | 0 | 1.0 | ambig_same_strand | spliced_annot | mRNA OSBPL7 ENST00000585051.1 |
| A00839:75:H7MFFDSXY:3:2130:5999:25567 | chr3:184,352,479-184,353,037 | mrna | 145M222N5M | -559 | 1 | 1 | 1.0 | ambig_same_strand | spliced_annot | mRNA CLCN2 ENST00000491162.1 |
| A00839:75:H7MFFDSXY:3:1108:8576:9204 | chr1:154,601,859-154,602,314 | mrna | 150M | -456 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA ADAR ENST00000680270.2 |
| A00839:75:H7MFFDSXY:3:2106:27145:12023 | chrX:109,375,769-109,376,246 | nrna | 150M | -478 | 1 | 0 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chrX:109,372,905-109,482,086(-) |
| A00839:75:H7MFFDSXY:3:2214:4761:31031 | chr10:102,869,391-102,869,697 | nrna | 119M36N31M | 307 | 1 | 0 | 1.0 | ambig_same_strand | spliced_unannot | nRNA chr10:102,854,258-102,901,899(+) |
| A00839:75:H7MFFDSXY:3:1669:18629:4225 | chr16:20,372,722-20,373,182 | nrna | 150M | -461 | 1 | 0 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chr16:20,359,174-20,404,737(-) |
| A00839:75:H7MFFDSXY:3:2627:2040:11600 | chr11:65,617,181-65,617,638 | nrna | 150M | 458 | 1 | 0 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chr11:65,615,775-65,637,439(+) |
| A00839:75:H7MFFDSXY:3:1561:32470:6151 | chr2:85,306,158-85,306,633 | nrna | 150M | 476 | 1 | 0 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chr2:85,133,391-85,310,387(+) |
| A00839:75:H7MFFDSXY:3:2358:28456:2519 | chr3:52,451,133-52,451,579 | nrna | 150M | -447 | 1 | 0 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chr3:52,451,099-52,454,041(-) |
| A00839:75:H7MFFDSXY:3:2418:20880:1454 | chr9:134,317,379-134,317,685 | nrna | 136M26N14M | 307 | 1 | 0 | 1.0 | ambig_same_strand | spliced_unannot | nRNA chr9:134,317,097-134,406,394(+) |
| A00839:75:H7MFFDSXY:3:1627:12292:24940 | chr7:128,837,346-128,837,784 | nrna | 150M | -439 | 1 | 0 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chr7:128,830,405-128,859,274(+) |
| A00839:75:H7MFFDSXY:3:1405:4643:18239 | chrX:41,342,042-41,342,822 | mrna | 150M | -781 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA DDX3X ENST00000615742.4 |

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

- Error categories: `results/vcap_no_mm_hotspots_2026-05-17/gdna_false_rna_categories.tsv`
- Hotspot windows: `results/vcap_no_mm_hotspots_2026-05-17/gdna_false_rna_windows.tsv`
- Hotspot EM loci: `results/vcap_no_mm_hotspots_2026-05-17/gdna_false_rna_loci.tsv`
- Assigned RNA targets: `results/vcap_no_mm_hotspots_2026-05-17/gdna_false_rna_targets.tsv`
- Representative fragments: `results/vcap_no_mm_hotspots_2026-05-17/sample_false_rna_fragments.tsv`
- High-confidence examples: `results/vcap_no_mm_hotspots_2026-05-17/high_confidence_false_rna_fragments.tsv`
- Summary JSON: `results/vcap_no_mm_hotspots_2026-05-17/summary.json`
