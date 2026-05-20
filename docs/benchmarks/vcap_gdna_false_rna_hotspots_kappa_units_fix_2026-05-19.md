# VCaP gDNA False-RNA Hotspot Analysis - 2026-05-16

BAM: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/kappa_units_fix/annotated.bam`

Truth source is derived from query-name flowcell ID. This analysis focuses on gDNA-source fragments (`H7MFFDSXY`) that Rigel assigned to an RNA pool (`mRNA` or `nRNA`). Counting rule is one primary read1 record per fragment.

Fragments counted: 32,218,601

gDNA-source fragments: 18,228,677

gDNA -> RNA false positives: 2,393,558 (13.13% of gDNA-source fragments)

False-positive split: 1,062,195 mRNA (44.38% of false RNA) and 1,331,363 nRNA (55.62% of false RNA).

## What Dominates

The false-positive RNA problem is concentrated in unspliced, genic ambiguous fragments. Across false RNA calls, 95.95% are `ZS=unspliced` and 96.94% are either `ambig_same_strand` or `ambig_opp_strand`.

The top 10 EM loci account for 13.38% of all gDNA -> RNA false positives. The top 10 genomic 10,000 bp windows account for 1.13%. So this is not only a few pathological coordinates; it is a recurrent ambiguous-region failure mode, with some strong hotspots.

## Error Categories

| Pred detail | ZC | ZS | Count | False RNA frac | gDNA source frac |
| --- | --- | --- | --- | --- | --- |
| nrna | ambig_same_strand | unspliced | 877,209 | 36.65% | 4.81% |
| mrna | ambig_same_strand | unspliced | 830,287 | 34.69% | 4.55% |
| nrna | ambig_opp_strand | unspliced | 342,908 | 14.33% | 1.88% |
| mrna | ambig_opp_strand | unspliced | 197,129 | 8.24% | 1.08% |
| nrna | ambig_same_strand | spliced_implicit | 42,010 | 1.76% | 0.23% |
| nrna | unambig | unspliced | 28,111 | 1.17% | 0.15% |
| nrna | multimapper | unspliced | 15,016 | 0.63% | 0.08% |
| nrna | multimapper | spliced_annot | 10,799 | 0.45% | 0.06% |
| mrna | ambig_same_strand | spliced_annot | 9,159 | 0.38% | 0.05% |
| nrna | ambig_opp_strand | spliced_implicit | 6,821 | 0.28% | 0.04% |
| mrna | ambig_same_strand | spliced_implicit | 6,509 | 0.27% | 0.04% |
| mrna | multimapper | unspliced | 5,931 | 0.25% | 0.03% |
| mrna | multimapper | spliced_annot | 4,643 | 0.19% | 0.03% |
| nrna | ambig_same_strand | spliced_unannot | 4,135 | 0.17% | 0.02% |
| nrna | multimapper | spliced_unannot | 3,643 | 0.15% | 0.02% |
| mrna | unambig | spliced_annot | 3,506 | 0.15% | 0.02% |
| mrna | ambig_same_strand | spliced_unannot | 1,822 | 0.08% | 0.01% |
| mrna | ambig_opp_strand | spliced_implicit | 1,669 | 0.07% | 0.01% |
| mrna | multimapper | spliced_unannot | 1,392 | 0.06% | 0.01% |
| nrna | ambig_opp_strand | spliced_unannot | 558 | 0.02% | 0.00% |
| nrna | unambig | spliced_unannot | 153 | 0.01% | 0.00% |
| mrna | ambig_opp_strand | spliced_unannot | 148 | 0.01% | 0.00% |

## Top 10kb Genomic Windows

| Region | False RNA | False rate | gDNA total | mRNA | nRNA | Top target | Top category |
| --- | --- | --- | --- | --- | --- | --- | --- |
| chr11:65,500,001-65,510,000 | 3,698 | 87.92% | 4,206 | 1,993 | 1,705 | nRNA ENST00000698129.1 | mrna/ambig_opp_strand/unspliced |
| chr1:152,350,001-152,360,000 | 3,117 | 67.61% | 4,610 | 1,449 | 1,668 | mRNA FLG2 ENST00000388718.5 | nrna/ambig_opp_strand/unspliced |
| chrX:67,540,001-67,550,000 | 2,875 | 43.66% | 6,585 | 2,710 | 165 | mRNA AR ENST00000396043.4 | mrna/ambig_same_strand/unspliced |
| chr11:64,220,001-64,230,000 | 2,827 | 48.67% | 5,809 | 760 | 2,067 | nRNA chr11:64,223,798-64,226,158(-) | nrna/ambig_same_strand/unspliced |
| chr2:178,560,001-178,570,000 | 2,665 | 67.49% | 3,949 | 179 | 2,486 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr2:178,570,001-178,580,000 | 2,481 | 66.75% | 3,717 | 166 | 2,315 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr2:178,530,001-178,540,000 | 2,447 | 66.71% | 3,668 | 305 | 2,142 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr11:64,230,001-64,240,000 | 2,445 | 52.96% | 4,617 | 182 | 2,263 | nRNA chr11:64,229,213-64,234,292(-) | nrna/ambig_opp_strand/unspliced |
| chr9:106,920,001-106,930,000 | 2,357 | 49.12% | 4,798 | 1,422 | 935 | mRNA ZNF462 ENST00000277225.10 | mrna/ambig_same_strand/unspliced |
| chr2:178,540,001-178,550,000 | 2,168 | 68.87% | 3,148 | 84 | 2,084 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr21:33,550,001-33,560,000 | 1,803 | 43.74% | 4,122 | 1,723 | 80 | mRNA SON ENST00000300278.8 | mrna/ambig_same_strand/unspliced |
| chr2:178,610,001-178,620,000 | 1,752 | 65.96% | 2,656 | 8 | 1,744 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr9:32,630,001-32,640,000 | 1,745 | 50.27% | 3,471 | 240 | 1,505 | nRNA ENST00000242310.4 | nrna/unambig/unspliced |
| chr9:76,700,001-76,710,000 | 1,739 | 78.37% | 2,219 | 962 | 777 | mRNA PRUNE2 ENST00000426088.5 | mrna/ambig_opp_strand/unspliced |
| chr2:178,590,001-178,600,000 | 1,641 | 65.22% | 2,516 | 58 | 1,583 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr2:178,550,001-178,560,000 | 1,615 | 67.40% | 2,396 | 54 | 1,561 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr11:64,110,001-64,120,000 | 1,592 | 39.19% | 4,062 | 1,013 | 579 | mRNA FLRT1 ENST00000246841.3 | mrna/ambig_opp_strand/unspliced |
| chr11:64,260,001-64,270,000 | 1,578 | 29.29% | 5,387 | 60 | 1,518 | nRNA chr11:64,251,529-64,267,923(+) | nrna/ambig_same_strand/unspliced |
| chr17:28,610,001-28,620,000 | 1,570 | 58.60% | 2,679 | 37 | 1,533 | nRNA chr17:28,611,931-28,617,377(+) | nrna/ambig_opp_strand/unspliced |
| chr11:66,330,001-66,340,000 | 1,512 | 47.15% | 3,207 | 82 | 1,430 | nRNA chr11:66,334,493-66,339,875(+) | nrna/ambig_opp_strand/unspliced |
| chr11:66,850,001-66,860,000 | 1,509 | 42.20% | 3,576 | 645 | 864 | nRNA chr11:66,857,882-66,860,475(-) | mrna/ambig_opp_strand/unspliced |
| chr1:152,210,001-152,220,000 | 1,462 | 66.61% | 2,195 | 685 | 777 | nRNA chr1:152,168,124-152,332,686(+) | nrna/ambig_opp_strand/unspliced |
| chr2:178,580,001-178,590,000 | 1,442 | 63.00% | 2,289 | 37 | 1,405 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr2:178,600,001-178,610,000 | 1,422 | 62.23% | 2,285 | 16 | 1,406 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr15:89,620,001-89,630,000 | 1,410 | 65.40% | 2,156 | 581 | 829 | mRNA TICRR ENST00000268138.12 | nrna/ambig_opp_strand/unspliced |
| chr2:185,790,001-185,800,000 | 1,386 | 65.97% | 2,101 | 679 | 707 | mRNA FSIP2 ENST00000424728.6 | nrna/ambig_opp_strand/unspliced |
| chr2:178,770,001-178,780,000 | 1,386 | 62.40% | 2,221 | 26 | 1,360 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr17:28,630,001-28,640,000 | 1,352 | 49.36% | 2,739 | 241 | 1,111 | nRNA chr17:28,614,881-28,645,454(-) | nrna/ambig_opp_strand/unspliced |
| chr11:64,280,001-64,290,000 | 1,320 | 29.60% | 4,459 | 529 | 791 | nRNA chr11:64,286,419-64,289,494(+) | nrna/ambig_same_strand/unspliced |
| chr11:64,250,001-64,260,000 | 1,307 | 28.69% | 4,556 | 16 | 1,291 | nRNA chr11:64,251,529-64,267,923(+) | nrna/ambig_same_strand/unspliced |
| chr2:178,720,001-178,730,000 | 1,296 | 56.57% | 2,291 | 13 | 1,283 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr2:178,740,001-178,750,000 | 1,280 | 59.59% | 2,148 | 29 | 1,251 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr13:109,780,001-109,790,000 | 1,270 | 32.81% | 3,871 | 1,260 | 10 | mRNA IRS2 ENST00000375856.5 | mrna/ambig_same_strand/unspliced |
| chr17:29,570,001-29,580,000 | 1,247 | 42.56% | 2,930 | 153 | 1,094 | nRNA chr17:29,560,546-29,707,090(+) | nrna/ambig_opp_strand/unspliced |
| chr2:178,710,001-178,720,000 | 1,246 | 57.71% | 2,159 | 22 | 1,224 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr2:178,730,001-178,740,000 | 1,234 | 57.26% | 2,155 | 28 | 1,206 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr13:23,330,001-23,340,000 | 1,208 | 35.39% | 3,413 | 1,208 | 0 | mRNA SACS ENST00000402364.1 | mrna/ambig_same_strand/unspliced |
| chr8:76,850,001-76,860,000 | 1,199 | 41.06% | 2,920 | 1,136 | 63 | mRNA ZFHX4 ENST00000518282.5 | mrna/ambig_same_strand/unspliced |
| chr17:28,570,001-28,580,000 | 1,197 | 43.48% | 2,753 | 104 | 1,093 | nRNA chr17:28,573,114-28,576,895(-) | nrna/ambig_same_strand/unspliced |
| chr4:87,610,001-87,620,000 | 1,179 | 52.35% | 2,252 | 626 | 553 | mRNA DSPP ENST00000651931.1 | mrna/ambig_opp_strand/unspliced |

## Local RNA Evidence Stratification

False RNA calls are not confined to windows with abundant RNA-source reads. This matters
because it separates two mechanisms: real local RNA expression can pull ambiguous gDNA
into RNA, but ambiguous gDNA can also self-seed RNA components in windows with little or
no RNA evidence.

| RNA-source fragments in 10kb window | Windows | False RNA | False RNA frac | Weighted false rate |
| --- | --- | --- | --- | --- |
| 0 | 22,391 | 303,089 | 12.66% | 11.07% |
| 1-100 | 30,962 | 704,475 | 29.43% | 11.97% |
| 101-1,000 | 16,895 | 881,390 | 36.82% | 15.45% |
| >1,000 | 2,975 | 504,604 | 21.08% | 21.57% |

## Top EM Loci

The region span is the observed fragment span among false-positive fragments in the locus, not the full locus extent. Components such as `chr2:3`, `chr1:3`, and `chr7:3` are chromosome-scale mega-components rather than precise biological loci, so the 10kb windows and target tables above are the sharper coordinates for assay triage.

| Locus | Observed region | False RNA | False frac | mRNA | nRNA | Same | Opp | Mean ZW | ZW >= .90 | Top target |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr2:3 | chr2:3,130,406-241,979,860 | 85,468 | 3.57% | 20,570 | 64,898 | 40,934 | 43,394 | 0.491 | 8.52% | nRNA chr2:178,525,988-178,804,642(-) |
| chr1:3 | chr1:131,229-248,593,775 | 39,107 | 1.63% | 13,163 | 25,944 | 24,931 | 11,333 | 0.521 | 8.56% | nRNA chr1:152,168,124-152,332,686(+) |
| chr7:3 | chr7:997,379-158,587,769 | 31,013 | 1.30% | 8,645 | 22,368 | 23,329 | 7,242 | 0.447 | 4.44% | nRNA chr7:141,995,878-142,106,747(+) |
| chr11:3 | chr11:124,137-130,295,870 | 29,275 | 1.22% | 9,518 | 19,757 | 15,955 | 13,003 | 0.531 | 13.15% | nRNA ENST00000698129.1 |
| chr3:3 | chr3:1,092,672-198,225,656 | 25,236 | 1.05% | 8,326 | 16,910 | 18,341 | 6,671 | 0.480 | 5.38% | nRNA chr3:38,039,204-38,124,025(+) |
| chr5:3 | chr5:5,422,637-181,331,826 | 24,819 | 1.04% | 7,369 | 17,450 | 13,109 | 10,865 | 0.549 | 12.26% | nRNA chr5:141,136,682-141,241,744(-) |
| chr9:3 | chr9:14,524-138,331,545 | 24,259 | 1.01% | 8,006 | 16,253 | 16,294 | 6,797 | 0.506 | 8.34% | nRNA chr9:110,365,247-110,579,741(-) |
| chr8:3 | chr8:120,387-140,451,471 | 23,069 | 0.96% | 6,116 | 16,953 | 16,935 | 5,948 | 0.482 | 5.75% | nRNA chr8:109,362,460-109,537,207(+) |
| chr4:3 | chr4:1,900,622-190,195,896 | 19,071 | 0.80% | 6,591 | 12,480 | 12,388 | 6,471 | 0.475 | 5.48% | nRNA chr4:78,057,322-78,544,269(+) |
| chr12:3 | chr12:15,172-133,204,855 | 18,931 | 0.79% | 3,999 | 14,932 | 15,064 | 3,589 | 0.451 | 2.94% | nRNA chr12:26,335,351-26,833,194(-) |
| chr17:3 | chr17:60,222-83,221,879 | 17,992 | 0.75% | 3,593 | 14,399 | 9,397 | 8,169 | 0.519 | 5.80% | nRNA chr17:28,614,881-28,645,454(-) |
| chr6:3 | chr6:95,134-170,745,948 | 15,190 | 0.63% | 4,289 | 10,901 | 11,911 | 3,001 | 0.433 | 3.31% | nRNA chr6:56,457,986-56,642,996(-) |
| chr15:3 | chr15:20,144,097-101,876,410 | 14,754 | 0.62% | 3,565 | 11,189 | 7,685 | 6,940 | 0.508 | 3.55% | nRNA chr15:81,332,346-81,362,937(+) |
| chr19:3 | chr19:156,842-58,599,378 | 8,871 | 0.37% | 3,871 | 5,000 | 3,643 | 5,110 | 0.528 | 7.85% | mRNA ZNF615 ENST00000618487.4 |
| chr7:883 | chr7:42,757,557-156,969,899 | 8,395 | 0.35% | 2,306 | 6,089 | 4,151 | 3,774 | 0.495 | 5.94% | nRNA chr7:100,197,473-100,214,387(+) |
| chr5:1380 | chr5:62,306,180-150,396,827 | 6,485 | 0.27% | 819 | 5,666 | 6,135 | 293 | 0.687 | 18.13% | nRNA ENST00000611950.1 |
| chr10:3 | chr10:5,765,702-133,468,290 | 6,375 | 0.27% | 2,177 | 4,198 | 5,805 | 384 | 0.514 | 9.32% | mRNA ZNF518A ENST00000624776.4 |
| chr17:867 | chr17:29,562,309-30,185,893 | 5,420 | 0.23% | 1,384 | 4,036 | 1,489 | 3,918 | 0.489 | 7.51% | nRNA chr17:29,560,546-29,707,090(+) |
| chr16:3 | chr16:1,510,650-90,175,600 | 5,134 | 0.21% | 2,579 | 2,555 | 3,189 | 1,601 | 0.509 | 17.96% | mRNA ZFHX3 ENST00000268489.10 |
| chr19:2804 | chr19:36,054,750-55,456,609 | 4,495 | 0.19% | 1,278 | 3,217 | 3,912 | 445 | 0.441 | 5.63% | nRNA chr19:42,325,634-42,378,769(+) |
| chr11:3225 | chr11:5,254,207-64,321,809 | 4,434 | 0.19% | 1,104 | 3,330 | 2,500 | 1,912 | 0.691 | 30.42% | nRNA chr11:64,318,120-64,321,811(+) |
| chr13:5941 | chr13:96,453,381-110,512,536 | 4,414 | 0.18% | 236 | 4,178 | 2,831 | 1,506 | 0.416 | 3.92% | nRNA chr13:110,307,283-110,504,595(+) |
| chrX:21523 | chrX:67,521,831-67,853,372 | 4,233 | 0.18% | 3,289 | 944 | 4,199 | 1 | 0.408 | 8.46% | mRNA AR ENST00000396043.4 |
| chr18:3 | chr18:1,272,297-65,880,930 | 3,902 | 0.16% | 1,302 | 2,600 | 2,075 | 1,784 | 0.494 | 13.22% | nRNA chr18:58,535,414-58,538,552(+) |
| chr14:3 | chr14:16,373,182-83,170,058 | 3,856 | 0.16% | 1,919 | 1,937 | 2,111 | 1,481 | 0.550 | 15.82% | mRNA ADAM20 ENST00000256389.5 |
| chr3:1097 | chr3:49,940,348-190,656,853 | 3,727 | 0.16% | 687 | 3,040 | 3,316 | 392 | 0.463 | 3.84% | nRNA chr3:52,495,337-52,524,495(+) |
| chr11:3926 | chr11:67,351,755-85,758,001 | 3,510 | 0.15% | 1,100 | 2,410 | 2,041 | 1,449 | 0.458 | 2.59% | mRNA SYTL2 ENST00000634661.1 |
| chr17:8543 | chr17:10,300,771-46,371,505 | 3,332 | 0.14% | 111 | 3,221 | 16 | 3,292 | 0.401 | 3.84% | nRNA chr17:10,383,155-10,537,862(+) |
| chrX:3 | chrX:7,049,962-155,568,141 | 3,308 | 0.14% | 1,348 | 1,960 | 2,181 | 1,015 | 0.508 | 6.50% | mRNA FRMPD4 ENST00000657982.1 |
| chr20:3 | chr20:13,849,252-64,292,591 | 3,267 | 0.14% | 1,821 | 1,446 | 2,348 | 864 | 0.622 | 27.27% | nRNA chr20:61,599,760-61,937,249(+) |
| chr21:13345 | chr21:11,628,714-33,931,473 | 3,256 | 0.14% | 2,288 | 968 | 2,532 | 694 | 0.429 | 3.10% | mRNA SON ENST00000300278.8 |
| chr3:14073 | chr3:9,750,038-195,443,031 | 3,121 | 0.13% | 1,447 | 1,674 | 1,482 | 1,616 | 0.486 | 3.14% | mRNA IGSF10 ENST00000282466.4 |
| chr11:3723 | chr11:59,712,933-64,328,765 | 2,946 | 0.12% | 72 | 2,874 | 2,869 | 40 | 0.735 | 9.27% | nRNA chr11:64,251,529-64,267,923(+) |
| chr11:3821 | chr11:63,974,620-64,166,101 | 2,790 | 0.12% | 1,324 | 1,466 | 927 | 1,857 | 0.483 | 1.72% | mRNA FLRT1 ENST00000246841.3 |
| chr21:4 | chr21:8,202,827-45,024,528 | 2,655 | 0.11% | 531 | 2,124 | 980 | 355 | 0.680 | 45.88% | nRNA chr21:8,380,664-8,450,668(+) |
| chr6:17233 | chr6:26,026,904-166,439,391 | 2,645 | 0.11% | 884 | 1,761 | 1,884 | 711 | 0.588 | 26.24% | nRNA chr6:26,457,949-26,469,653(+) |
| chr17:247 | chr17:15,995,681-81,720,975 | 2,576 | 0.11% | 889 | 1,687 | 2,069 | 462 | 0.398 | 4.00% | nRNA chr17:81,395,456-81,466,331(+) |
| chr12:4494 | chr12:1,946,055-96,952,309 | 2,553 | 0.11% | 544 | 2,009 | 1,313 | 1,219 | 0.462 | 1.61% | nRNA chr12:2,820,087-2,859,791(+) |
| chr15:7093 | chr15:41,183,673-99,133,452 | 2,339 | 0.10% | 1,047 | 1,292 | 1,221 | 1,109 | 0.624 | 10.13% | mRNA SYNM ENST00000336292.11 |
| chr8:19853 | chr8:76,703,819-92,818,668 | 2,281 | 0.10% | 2,036 | 245 | 2,279 | 1 | 0.656 | 36.52% | mRNA ZFHX4 ENST00000518282.5 |

## Top Assigned RNA Targets

| Target | Observed region | False RNA | False frac | mRNA | nRNA | Same | Opp | Mean ZW |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| nRNA chr2:178,525,988-178,804,642(-) | chr2:32,916,265-178,804,640 | 14,066 | 0.59% | 0 | 14,066 | 417 | 13,623 | 0.735 |
| nRNA chr2:178,523,417-178,776,534(+) | chr2:178,526,721-178,776,529 | 6,812 | 0.28% | 0 | 6,812 | 0 | 6,801 | 0.452 |
| nRNA chr11:64,251,529-64,267,923(+) | chr11:64,251,522-64,328,765 | 2,809 | 0.12% | 0 | 2,809 | 2,800 | 0 | 0.742 |
| nRNA chr9:110,365,247-110,579,741(-) | chr9:110,366,202-110,579,718 | 2,170 | 0.09% | 0 | 2,170 | 2,166 | 0 | 0.650 |
| nRNA chr17:29,560,546-29,707,090(+) | chr17:29,562,309-29,703,236 | 1,750 | 0.07% | 0 | 1,750 | 226 | 1,519 | 0.482 |
| nRNA chr21:8,380,664-8,450,668(+) | mixed:101,586-45,024,528 | 1,715 | 0.07% | 0 | 1,715 | 137 | 0 | 0.829 |
| mRNA ZFHX4 ENST00000518282.5 | chr8:76,704,045-92,818,668 | 1,711 | 0.07% | 1,711 | 0 | 1,710 | 1 | 0.767 |
| nRNA chr5:141,136,682-141,241,744(-) | chr5:141,137,019-141,241,485 | 1,704 | 0.07% | 0 | 1,704 | 3 | 1,677 | 0.640 |
| nRNA ENST00000698129.1 | chr11:65,499,302-65,506,876 | 1,701 | 0.07% | 0 | 1,701 | 0 | 1,692 | 0.846 |
| nRNA chr1:152,168,124-152,332,686(+) | chr1:152,179,714-152,326,058 | 1,688 | 0.07% | 0 | 1,688 | 13 | 1,246 | 0.796 |
| mRNA MUC16 ENST00000397910.8 | chr19:8,886,745-8,981,296 | 1,670 | 0.07% | 1,670 | 0 | 1,660 | 0 | 0.518 |
| mRNA FSIP2 ENST00000424728.6 | chr2:185,788,656-185,833,291 | 1,618 | 0.07% | 1,618 | 0 | 576 | 1,039 | 0.782 |
| nRNA chr2:178,528,739-178,633,078(+) | chr2:178,528,814-178,632,987 | 1,575 | 0.07% | 0 | 1,575 | 0 | 1,575 | 0.138 |
| nRNA chr8:109,362,460-109,537,207(+) | chr8:109,362,473-109,530,370 | 1,569 | 0.07% | 0 | 1,569 | 1,558 | 0 | 0.569 |
| nRNA chr7:141,995,878-142,106,747(+) | chr7:142,005,153-142,106,142 | 1,558 | 0.07% | 0 | 1,558 | 1,533 | 0 | 0.668 |
| mRNA AR ENST00000396043.4 | chrX:67,544,227-67,724,277 | 1,536 | 0.06% | 1,536 | 0 | 1,528 | 1 | 0.492 |
| mRNA FLG2 ENST00000388718.5 | chr1:152,350,091-152,357,632 | 1,449 | 0.06% | 1,449 | 0 | 0 | 1,433 | 0.750 |
| mRNA AHNAK ENST00000378024.9 | chr11:62,516,582-62,534,062 | 1,433 | 0.06% | 1,433 | 0 | 1,423 | 0 | 0.689 |
| nRNA chr1:152,335,378-152,364,080(+) | chr1:152,344,441-152,359,031 | 1,432 | 0.06% | 0 | 1,432 | 1 | 1,405 | 0.800 |
| nRNA ENST00000242310.4 | chr9:32,629,863-32,635,660 | 1,363 | 0.06% | 0 | 1,363 | 0 | 543 | 0.846 |
| nRNA chr6:33,621,321-33,696,574(+) | chr6:33,621,491-33,892,590 | 1,361 | 0.06% | 0 | 1,361 | 1,317 | 34 | 0.519 |
| nRNA chr11:64,823,051-64,844,653(-) | chr11:64,687,206-64,844,644 | 1,357 | 0.06% | 0 | 1,357 | 1,355 | 0 | 0.650 |
| nRNA chr2:151,485,335-151,734,487(-) | chr2:151,485,504-151,733,337 | 1,344 | 0.06% | 0 | 1,344 | 1,225 | 115 | 0.418 |
| mRNA SON ENST00000300278.8 | chr21:33,549,517-33,560,460 | 1,313 | 0.05% | 1,313 | 0 | 1,255 | 57 | 0.607 |
| nRNA chr17:28,614,881-28,645,454(-) | chr17:28,614,977-28,645,427 | 1,305 | 0.05% | 0 | 1,305 | 955 | 347 | 0.695 |
| mRNA FAT1 ENST00000441802.7 | chr4:83,052,074-186,709,837 | 1,304 | 0.05% | 1,304 | 0 | 1,300 | 1 | 0.796 |
| mRNA IRS2 ENST00000375856.5 | chr13:109,755,854-109,786,239 | 1,301 | 0.05% | 1,301 | 0 | 1,142 | 159 | 0.818 |
| nRNA chr11:64,229,213-64,234,292(-) | chr11:64,229,221-64,234,291 | 1,292 | 0.05% | 0 | 1,292 | 0 | 1,286 | 0.796 |
| nRNA chr17:28,455,751-28,614,197(-) | chr17:28,469,436-28,614,194 | 1,242 | 0.05% | 0 | 1,242 | 451 | 789 | 0.479 |
| mRNA BSN ENST00000296452.5 | chr3:49,624,972-49,671,352 | 1,227 | 0.05% | 1,227 | 0 | 1,223 | 0 | 0.767 |
| nRNA chr2:178,541,746-178,597,681(+) | chr2:178,542,130-178,597,679 | 1,221 | 0.05% | 0 | 1,221 | 0 | 1,220 | 0.168 |
| mRNA PCLO ENST00000423517.6 | chr7:82,821,985-83,162,900 | 1,205 | 0.05% | 1,205 | 0 | 1,197 | 0 | 0.617 |
| nRNA chr8:102,254,197-102,412,759(-) | chr8:102,254,219-102,412,422 | 1,198 | 0.05% | 0 | 1,198 | 1,179 | 13 | 0.555 |
| nRNA chr13:110,307,283-110,504,595(+) | chr13:110,307,666-110,504,548 | 1,168 | 0.05% | 0 | 1,168 | 845 | 287 | 0.309 |
| mRNA FAT4 ENST00000394329.9 | chr4:125,315,073-125,492,011 | 1,166 | 0.05% | 1,166 | 0 | 1,164 | 0 | 0.716 |
| mRNA SACS ENST00000402364.1 | chr13:23,329,055-23,375,270 | 1,163 | 0.05% | 1,163 | 0 | 1,163 | 0 | 0.610 |
| nRNA chr13:110,149,060-110,299,118(-) | chr13:110,149,066-110,298,805 | 1,124 | 0.05% | 0 | 1,124 | 1,113 | 0 | 0.309 |
| mRNA CCDC168 ENST00000322527.4 | chr13:102,729,396-102,759,046 | 1,113 | 0.05% | 1,113 | 0 | 1,113 | 0 | 0.341 |
| nRNA chr4:78,057,322-78,544,269(+) | chr4:78,057,885-78,541,268 | 1,093 | 0.05% | 0 | 1,093 | 1,088 | 0 | 0.364 |
| nRNA chr17:10,383,155-10,537,862(+) | chr17:10,390,185-10,537,841 | 1,058 | 0.04% | 0 | 1,058 | 1 | 1,051 | 0.407 |

## Representative False-Positive Fragments From Top Loci

| Locus | QNAME | Region | Detail | CIGAR | TLEN | NH | NM | ZW | ZC | ZS | Target |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr2:3 | A00839:75:H7MFFDSXY:3:2324:24279:31955 | chr2:202,765,716-202,765,969 | mrna | 150M | -254 | 1 | 0 | 0.9526009559631348 | ambig_same_strand | unspliced | mRNA FAM117B ENST00000392238.3 |
| chr2:3 | A00839:75:H7MFFDSXY:3:1619:12436:22435 | chr2:46,369,650-46,369,895 | nrna | 150M | -246 | 1 | 0 | 0.5015048980712891 | ambig_same_strand | unspliced | nRNA chr2:46,297,406-46,386,697(+) |
| chr2:3 | A00839:75:H7MFFDSXY:3:1478:11360:25175 | chr2:17,721,060-17,721,290 | nrna | 150M | 231 | 1 | 0 | 0.2741115689277649 | ambig_same_strand | unspliced | nRNA chr2:17,665,446-17,800,195(-) |
| chr2:3 | A00839:75:H7MFFDSXY:3:1552:32497:13463 | chr2:178,718,006-178,718,262 | nrna | 150M | -257 | 1 | 0 | 0.09118147194385529 | ambig_opp_strand | unspliced | nRNA chr2:178,528,764-178,738,038(+) |
| chr1:3 | A00839:75:H7MFFDSXY:3:2324:7111:31031 | chr1:156,595,150-156,595,562 | nrna | 150M | -413 | 2 | 0 | 0.0578426867723465 | multimapper | unspliced | nRNA chr1:156,591,755-156,595,741(+) |
| chr1:3 | A00839:75:H7MFFDSXY:3:1376:18304:15781 | chr1:109,750,676-109,750,856 | nrna | 150M | 181 | 1 | 0 | 0.639573335647583 | ambig_opp_strand | unspliced | nRNA chr1:109,750,355-109,751,546(-) |
| chr1:3 | A00839:75:H7MFFDSXY:3:1436:25536:12712 | chr1:152,350,515-152,350,871 | mrna | 150M | 357 | 1 | 0 | 0.5728872418403625 | ambig_opp_strand | unspliced | mRNA FLG2 ENST00000388718.5 |
| chr1:3 | A00839:75:H7MFFDSXY:3:2130:1289:21762 | chr1:152,353,387-152,353,563 | nrna | 150M | -177 | 1 | 1 | 0.9229825139045715 | ambig_opp_strand | unspliced | nRNA chr1:152,335,378-152,364,080(+) |
| chr7:3 | A00839:75:H7MFFDSXY:3:2571:8169:34898 | chr7:117,297,519-117,297,759 | nrna | 150M | 241 | 1 | 0 | 0.6906374096870422 | ambig_opp_strand | unspliced | nRNA chr7:117,277,530-117,323,114(-) |
| chr7:3 | A00839:75:H7MFFDSXY:3:2466:19416:14387 | chr7:154,872,565-154,872,920 | nrna | 150M | -356 | 1 | 2 | 0.053429849445819855 | ambig_same_strand | unspliced | nRNA chr7:154,868,055-154,876,386(+) |
| chr7:3 | A00839:75:H7MFFDSXY:3:1309:10086:6308 | chr7:108,232,308-108,232,482 | nrna | 150M | 175 | 1 | 1 | 0.6221303939819336 | ambig_same_strand | unspliced | nRNA chr7:108,148,116-108,456,339(-) |
| chr7:3 | A00839:75:H7MFFDSXY:3:1666:13919:21277 | chr7:90,163,010-90,163,350 | nrna | 150M | -341 | 1 | 0 | 0.380150705575943 | ambig_opp_strand | unspliced | nRNA chr7:90,161,265-90,164,829(+) |
| chr11:3 | A00839:75:H7MFFDSXY:3:2224:27778:12054 | chr11:65,223,531-65,223,708 | mrna | 150M | -178 | 1 | 0 | 0.5031686425209045 | ambig_same_strand | unspliced | mRNA ENSG00000293475 ENST00000534471.1 |
| chr11:3 | A00839:75:H7MFFDSXY:3:1234:31620:12477 | chr11:64,219,947-64,220,282 | nrna | 150M | -336 | 1 | 0 | 0.38957709074020386 | ambig_same_strand | unspliced | nRNA chr11:64,219,509-64,220,669(+) |
| chr11:3 | A00839:75:H7MFFDSXY:3:2675:25003:35055 | chr11:64,225,500-64,225,802 | nrna | 150M | 303 | 1 | 0 | 0.16965807974338531 | ambig_same_strand | unspliced | nRNA chr11:64,223,852-64,225,870(-) |
| chr11:3 | A00839:75:H7MFFDSXY:3:1624:11686:30154 | chr11:65,500,844-65,501,189 | mrna | 150M | -346 | 1 | 0 | 0.7094918489456177 | ambig_opp_strand | unspliced | mRNA MALAT1 ENST00000710933.1 |
| chr3:3 | A00839:75:H7MFFDSXY:3:2437:16893:33786 | chr3:47,343,515-47,343,751 | nrna | 150M | -237 | 1 | 0 | 0.6793757677078247 | ambig_same_strand | unspliced | nRNA chr3:47,282,943-47,344,577(+) |
| chr3:3 | A00839:75:H7MFFDSXY:3:2619:32958:33332 | chr3:53,775,974-53,776,112 | nrna | 139M | -139 | 1 | 1 | 0.1555773764848709 | ambig_same_strand | unspliced | nRNA chr3:53,494,610-53,781,924(+) |
| chr3:3 | A00839:75:H7MFFDSXY:3:1233:16206:23296 | chr3:53,800,253-53,800,507 | nrna | 150M | -255 | 1 | 0 | 0.5940913558006287 | ambig_same_strand | unspliced | nRNA chr3:53,799,890-53,813,712(+) |
| chr3:3 | A00839:75:H7MFFDSXY:3:2552:23764:25801 | chr3:47,891,475-47,891,772 | mrna | 150M | 298 | 1 | 1 | 0.7096536159515381 | ambig_same_strand | unspliced | mRNA MAP4 ENST00000497735.1 |
| chr5:3 | A00839:75:H7MFFDSXY:3:1163:14262:25441 | chr5:141,183,091-141,183,317 | nrna | 150M | 227 | 1 | 0 | 0.8352029323577881 | ambig_opp_strand | unspliced | nRNA chr5:141,136,682-141,241,744(-) |
| chr5:3 | A00839:75:H7MFFDSXY:3:1512:7401:2409 | chr5:141,235,450-141,235,769 | nrna | 150M | -320 | 1 | 0 | 0.588782548904419 | ambig_opp_strand | unspliced | nRNA ENST00000524813.2 |
| chr5:3 | A00839:75:H7MFFDSXY:3:2438:24144:12837 | chr5:140,652,894-140,653,167 | nrna | 150M | -274 | 1 | 0 | 0.5811043977737427 | ambig_same_strand | unspliced | nRNA chr5:140,647,828-140,662,480(+) |
| chr5:3 | A00839:75:H7MFFDSXY:3:2466:28718:20353 | chr5:115,391,605-115,391,769 | nrna | 150M | 165 | 1 | 0 | 0.697519838809967 | ambig_same_strand | unspliced | nRNA ENST00000515849.2 |
| chr9:3 | A00839:75:H7MFFDSXY:3:1647:25554:35289 | chr9:32,633,742-32,634,092 | nrna | 150M | 351 | 1 | 0 | 0.7445800304412842 | ambig_opp_strand | unspliced | nRNA ENST00000242310.4 |
| chr9:3 | A00839:75:H7MFFDSXY:3:2504:2248:3568 | chr9:92,615,947-92,616,184 | mrna | 150M | -238 | 1 | 0 | 0.805242121219635 | ambig_opp_strand | unspliced | mRNA CENPP ENST00000375587.8 |
| chr9:3 | A00839:75:H7MFFDSXY:3:2512:9552:36010 | chr9:106,924,052-106,924,272 | nrna | 150M | -221 | 1 | 1 | 0.9887856245040894 | ambig_same_strand | unspliced | nRNA chr9:106,863,165-106,930,575(+) |
| chr9:3 | A00839:75:H7MFFDSXY:3:2358:21423:33833 | chr9:110,481,103-110,481,314 | nrna | 150M | 212 | 1 | 0 | 0.872234046459198 | ambig_same_strand | unspliced | nRNA chr9:110,365,247-110,579,741(-) |
| chr8:3 | A00839:75:H7MFFDSXY:3:1517:4110:30279 | chr8:99,891,808-99,892,048 | nrna | 150M | 241 | 1 | 1 | 0.03282276913523674 | ambig_same_strand | unspliced | nRNA chr8:99,874,995-99,893,707(-) |
| chr8:3 | A00839:75:H7MFFDSXY:3:1343:12310:18959 | chr8:72,937,822-72,937,992 | mrna | 150M | -171 | 1 | 0 | 0.9480561017990112 | ambig_same_strand | unspliced | mRNA KCNB2 ENST00000523207.2 |
| chr8:3 | A00839:75:H7MFFDSXY:3:2468:24605:10848 | chr8:124,143,311-124,143,639 | nrna | 150M | 329 | 1 | 2 | 0.0920981615781784 | ambig_opp_strand | unspliced | nRNA chr8:124,058,536-124,171,451(-) |
| chr8:3 | A00839:75:H7MFFDSXY:3:1643:7166:6731 | chr8:70,660,679-70,660,975 | nrna | 150M | 297 | 1 | 0 | 0.39539653062820435 | ambig_opp_strand | unspliced | nRNA chr8:70,637,265-70,669,185(-) |
| chr4:3 | A00839:75:H7MFFDSXY:3:1168:30309:18067 | chr4:99,949,635-99,949,866 | nrna | 150M | 232 | 1 | 0 | 0.8801082968711853 | ambig_same_strand | unspliced | nRNA chr4:99,948,777-99,950,268(-) |
| chr4:3 | A00839:75:H7MFFDSXY:3:2346:21956:29465 | chr4:87,612,481-87,612,737 | mrna | 150M | -257 | 1 | 0 | 0.8727400302886963 | ambig_opp_strand | unspliced | mRNA DSPP ENST00000651931.1 |
| chr4:3 | A00839:75:H7MFFDSXY:3:1323:31060:13604 | chr4:147,879,041-147,879,416 | nrna | 150M | -376 | 1 | 0 | 0.013599333353340626 | ambig_same_strand | unspliced | nRNA chr4:147,860,982-148,072,776(+) |
| chr4:3 | A00839:75:H7MFFDSXY:3:2524:16685:17456 | chr4:72,288,662-72,288,942 | nrna | 150M | 281 | 1 | 0 | 0.6099640130996704 | ambig_same_strand | unspliced | nRNA chr4:72,280,968-72,569,221(-) |
| chr12:3 | A00839:75:H7MFFDSXY:3:1447:17409:25786 | chr12:39,643,784-39,644,034 | nrna | 150M | -251 | 1 | 1 | 0.6637797355651855 | ambig_same_strand | unspliced | nRNA chr12:39,626,166-39,721,914(+) |
| chr12:3 | A00839:75:H7MFFDSXY:3:2642:1316:35211 | chr12:39,646,372-39,646,580 | nrna | 150M | -209 | 1 | 0 | 0.6949716806411743 | ambig_same_strand | unspliced | nRNA chr12:39,626,166-39,721,914(+) |
| chr12:3 | A00839:75:H7MFFDSXY:3:2432:25988:15060 | chr12:51,379,677-51,391,603 | mrna | 147M11716N3M | -11927 | 2 | 0 | 0.9868096709251404 | multimapper | spliced_annot | mRNA GALNT6 ENST00000604381.5 |
| chr12:3 | A00839:75:H7MFFDSXY:3:2408:15799:15796 | chr12:21,476,964-21,477,166 | nrna | 96M1I53M | 203 | 1 | 2 | 0.7313284277915955 | ambig_same_strand | unspliced | nRNA chr12:21,469,801-21,501,551(-) |

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
| A00839:75:H7MFFDSXY:3:2501:10402:30029 | chrX:136,346,372-136,346,713 | nrna | 150M | -342 | 8 | 0 | 1.0 | multimapper | unspliced | nRNA chrX:104,566,198-105,767,829(+) |
| A00839:75:H7MFFDSXY:3:2316:1985:27571 | chr1:153,628,145-153,628,860 | mrna | 150M | -716 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA S100A1 ENST00000368698.3 |
| A00839:75:H7MFFDSXY:3:1616:24542:14184 | chr7:121,329,241-121,329,632 | mrna | 147M179N3M | 392 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA WNT16 ENST00000222462.3 |
| A00839:75:H7MFFDSXY:3:1262:16161:19993 | chr6:30,345,029-30,345,458 | nrna | 150M | 430 | 1 | 0 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chr6:30,328,888-30,346,857(+) |
| A00839:75:H7MFFDSXY:3:1134:9652:27665 | chrX:49,122,360-49,122,829 | nrna | 150M | -470 | 1 | 1 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chrX:49,113,406-49,123,735(-) |
| A00839:75:H7MFFDSXY:3:1317:17571:2675 | chr9:35,703,459-35,704,021 | nrna | 150M | 563 | 1 | 0 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chr9:35,696,947-35,732,195(-) |
| A00839:75:H7MFFDSXY:3:2369:27136:32675 | chr12:112,134,727-112,135,072 | mrna | 1S149M | -346 | 1 | 2 | 1.0 | ambig_same_strand | spliced_annot | mRNA TRAFD1 ENST00000548092.5 |
| A00839:75:H7MFFDSXY:3:1549:11324:7983 | chrX:20,130,595-20,138,599 | mrna | 89M2697N61M | -8005 | 1 | 0 | 1.0 | unambig | spliced_annot | mRNA EIF1AX ENST00000379607.10 |
| A00839:75:H7MFFDSXY:3:2119:18204:22842 | chr2:118,088,664-118,096,543 | mrna | 3M7752N100M1D24M | 7880 | 1 | 1 | 1.0 | unambig | spliced_annot | mRNA INSIG2 ENST00000467223.5 |
| A00839:75:H7MFFDSXY:3:2121:10420:34350 | chr3:49,795,593-49,796,148 | nrna | 150M | -556 | 1 | 0 | 1.0 | ambig_same_strand | spliced_implicit | nRNA chr3:49,790,731-49,799,873(-) |

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

- Error categories: `results/vcap_gdna_false_rna_hotspots_kappa_units_fix_2026-05-19/gdna_false_rna_categories.tsv`
- Hotspot windows: `results/vcap_gdna_false_rna_hotspots_kappa_units_fix_2026-05-19/gdna_false_rna_windows.tsv`
- Hotspot EM loci: `results/vcap_gdna_false_rna_hotspots_kappa_units_fix_2026-05-19/gdna_false_rna_loci.tsv`
- Assigned RNA targets: `results/vcap_gdna_false_rna_hotspots_kappa_units_fix_2026-05-19/gdna_false_rna_targets.tsv`
- Representative fragments: `results/vcap_gdna_false_rna_hotspots_kappa_units_fix_2026-05-19/sample_false_rna_fragments.tsv`
- High-confidence examples: `results/vcap_gdna_false_rna_hotspots_kappa_units_fix_2026-05-19/high_confidence_false_rna_fragments.tsv`
- Summary JSON: `results/vcap_gdna_false_rna_hotspots_kappa_units_fix_2026-05-19/summary.json`
