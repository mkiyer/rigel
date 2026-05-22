# VCaP gDNA False-RNA Hotspot Analysis - 2026-05-20

BAM: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/exon_strand_deconv_v1/annotated.bam`

Truth source is derived from query-name flowcell ID. This analysis focuses on gDNA-source fragments (`H7MFFDSXY`) that Rigel assigned to an RNA pool (`mRNA` or `nRNA`). Counting rule is one primary read1 record per fragment.

Fragments counted: 32,218,601

gDNA-source fragments: 18,228,677

gDNA -> RNA false positives: 2,500,290 (13.72% of gDNA-source fragments)

False-positive split: 1,171,473 mRNA (46.85% of false RNA) and 1,328,817 nRNA (53.15% of false RNA).

## What Dominates

The false-positive RNA problem is concentrated in unspliced, genic ambiguous fragments. Across false RNA calls, 96.12% are `ZS=unspliced` and 96.97% are either `ambig_same_strand` or `ambig_opp_strand`.

The top 10 EM loci account for 13.20% of all gDNA -> RNA false positives. The top 10 genomic 10,000 bp windows account for 1.13%. So this is not only a few pathological coordinates; it is a recurrent ambiguous-region failure mode, with some strong hotspots.

## Error Categories

| Pred detail | ZC | ZS | Count | False RNA frac | gDNA source frac |
| --- | --- | --- | --- | --- | --- |
| mrna | ambig_same_strand | unspliced | 918,169 | 36.72% | 5.04% |
| nrna | ambig_same_strand | unspliced | 854,615 | 34.18% | 4.69% |
| nrna | ambig_opp_strand | unspliced | 361,080 | 14.44% | 1.98% |
| mrna | ambig_opp_strand | unspliced | 217,716 | 8.71% | 1.19% |
| nrna | ambig_same_strand | spliced_implicit | 41,919 | 1.68% | 0.23% |
| nrna | unambig | unspliced | 30,380 | 1.22% | 0.17% |
| nrna | multimapper | unspliced | 14,926 | 0.60% | 0.08% |
| nrna | multimapper | spliced_annot | 10,662 | 0.43% | 0.06% |
| mrna | ambig_same_strand | spliced_annot | 9,159 | 0.37% | 0.05% |
| nrna | ambig_opp_strand | spliced_implicit | 6,821 | 0.27% | 0.04% |
| mrna | ambig_same_strand | spliced_implicit | 6,600 | 0.26% | 0.04% |
| mrna | multimapper | unspliced | 6,437 | 0.26% | 0.04% |
| mrna | multimapper | spliced_annot | 4,780 | 0.19% | 0.03% |
| nrna | ambig_same_strand | spliced_unannot | 4,108 | 0.16% | 0.02% |
| nrna | multimapper | spliced_unannot | 3,598 | 0.14% | 0.02% |
| mrna | unambig | spliced_annot | 3,506 | 0.14% | 0.02% |
| mrna | ambig_same_strand | spliced_unannot | 1,849 | 0.07% | 0.01% |
| mrna | ambig_opp_strand | spliced_implicit | 1,669 | 0.07% | 0.01% |
| mrna | multimapper | spliced_unannot | 1,437 | 0.06% | 0.01% |
| nrna | ambig_opp_strand | spliced_unannot | 555 | 0.02% | 0.00% |
| nrna | unambig | spliced_unannot | 153 | 0.01% | 0.00% |
| mrna | ambig_opp_strand | spliced_unannot | 151 | 0.01% | 0.00% |

## Top 10kb Genomic Windows

| Region | False RNA | False rate | gDNA total | mRNA | nRNA | Top target | Top category |
| --- | --- | --- | --- | --- | --- | --- | --- |
| chr11:65,500,001-65,510,000 | 3,809 | 90.56% | 4,206 | 1,990 | 1,819 | nRNA ENST00000698129.1 | mrna/ambig_opp_strand/unspliced |
| chr1:152,350,001-152,360,000 | 3,297 | 71.52% | 4,610 | 1,581 | 1,716 | mRNA FLG2 ENST00000388718.5 | nrna/ambig_opp_strand/unspliced |
| chr11:64,220,001-64,230,000 | 2,984 | 51.37% | 5,809 | 781 | 2,203 | nRNA chr11:64,218,896-64,223,857(+) | nrna/ambig_same_strand/unspliced |
| chrX:67,540,001-67,550,000 | 2,890 | 43.89% | 6,585 | 2,797 | 93 | mRNA AR ENST00000396043.4 | mrna/ambig_same_strand/unspliced |
| chr2:178,560,001-178,570,000 | 2,833 | 71.74% | 3,949 | 232 | 2,601 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr2:178,570,001-178,580,000 | 2,665 | 71.70% | 3,717 | 162 | 2,503 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr2:178,530,001-178,540,000 | 2,604 | 70.99% | 3,668 | 352 | 2,252 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr11:64,230,001-64,240,000 | 2,527 | 54.73% | 4,617 | 180 | 2,347 | nRNA chr11:64,229,213-64,234,292(-) | nrna/ambig_opp_strand/unspliced |
| chr9:106,920,001-106,930,000 | 2,411 | 50.25% | 4,798 | 2,036 | 375 | mRNA ZNF462 ENST00000277225.10 | mrna/ambig_same_strand/unspliced |
| chr2:178,540,001-178,550,000 | 2,282 | 72.49% | 3,148 | 112 | 2,170 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr2:178,610,001-178,620,000 | 1,880 | 70.78% | 2,656 | 15 | 1,865 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr21:33,550,001-33,560,000 | 1,847 | 44.81% | 4,122 | 1,784 | 63 | mRNA SON ENST00000300278.8 | mrna/ambig_same_strand/unspliced |
| chr9:32,630,001-32,640,000 | 1,829 | 52.69% | 3,471 | 289 | 1,540 | nRNA ENST00000242310.4 | nrna/unambig/unspliced |
| chr9:76,700,001-76,710,000 | 1,752 | 78.95% | 2,219 | 989 | 763 | nRNA chr9:76,704,847-76,710,661(+) | mrna/ambig_opp_strand/unspliced |
| chr2:178,590,001-178,600,000 | 1,735 | 68.96% | 2,516 | 47 | 1,688 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr11:64,110,001-64,120,000 | 1,705 | 41.97% | 4,062 | 1,160 | 545 | mRNA FLRT1 ENST00000246841.3 | mrna/ambig_opp_strand/unspliced |
| chr2:178,550,001-178,560,000 | 1,676 | 69.95% | 2,396 | 66 | 1,610 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr11:64,260,001-64,270,000 | 1,645 | 30.54% | 5,387 | 60 | 1,585 | nRNA chr11:64,251,529-64,267,923(+) | nrna/ambig_same_strand/unspliced |
| chr11:66,850,001-66,860,000 | 1,642 | 45.92% | 3,576 | 648 | 994 | nRNA chr11:66,857,882-66,860,475(-) | nrna/ambig_opp_strand/unspliced |
| chr11:66,330,001-66,340,000 | 1,623 | 50.61% | 3,207 | 78 | 1,545 | nRNA chr11:66,334,493-66,339,875(+) | nrna/ambig_opp_strand/unspliced |
| chr17:28,610,001-28,620,000 | 1,615 | 60.28% | 2,679 | 38 | 1,577 | nRNA chr17:28,611,931-28,617,377(+) | nrna/ambig_opp_strand/unspliced |
| chr2:178,600,001-178,610,000 | 1,565 | 68.49% | 2,285 | 27 | 1,538 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr1:152,210,001-152,220,000 | 1,552 | 70.71% | 2,195 | 755 | 797 | mRNA HRNR ENST00000368801.4 | nrna/ambig_opp_strand/unspliced |
| chr2:178,580,001-178,590,000 | 1,546 | 67.54% | 2,289 | 37 | 1,509 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr2:178,770,001-178,780,000 | 1,511 | 68.03% | 2,221 | 48 | 1,463 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr2:185,790,001-185,800,000 | 1,488 | 70.82% | 2,101 | 752 | 736 | mRNA FSIP2 ENST00000424728.6 | mrna/ambig_opp_strand/unspliced |
| chr15:89,620,001-89,630,000 | 1,483 | 68.78% | 2,156 | 569 | 914 | mRNA TICRR ENST00000268138.12 | nrna/ambig_opp_strand/unspliced |
| chr2:178,720,001-178,730,000 | 1,438 | 62.77% | 2,291 | 14 | 1,424 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |
| chr17:28,630,001-28,640,000 | 1,438 | 52.50% | 2,739 | 239 | 1,199 | nRNA chr17:28,614,881-28,645,454(-) | nrna/ambig_opp_strand/unspliced |
| chr2:178,730,001-178,740,000 | 1,367 | 63.43% | 2,155 | 33 | 1,334 | nRNA chr2:178,525,988-178,804,642(-) | nrna/ambig_opp_strand/unspliced |

## Local RNA Evidence Stratification

False RNA calls are not confined to windows with abundant RNA-source reads. This matters
because it separates two mechanisms: real local RNA expression can pull ambiguous gDNA
into RNA, but ambiguous gDNA can also self-seed RNA components in windows with little or
no RNA evidence.

| RNA-source fragments in 10kb window | Windows | False RNA | False RNA frac | Weighted false rate |
| --- | --- | --- | --- | --- |
| 0 | 20,448 | 315,027 | 12.60% | 11.66% |
| 1-100 | 29,954 | 733,356 | 29.33% | 12.60% |
| 101-1,000 | 16,746 | 922,118 | 36.88% | 16.22% |
| >1,000 | 2,976 | 529,789 | 21.19% | 22.65% |

## Top EM Loci

The region span is the observed fragment span among false-positive fragments in the locus, not the full locus extent. Components such as `chr2:3`, `chr1:3`, and `chr7:3` are chromosome-scale mega-components rather than precise biological loci, so the 10kb windows and target tables above are the sharper coordinates for assay triage.

| Locus | Observed region | False RNA | False frac | mRNA | nRNA | Same | Opp | Mean ZW | ZW >= .90 | Top target |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr2:3 | chr2:3,130,406-241,979,860 | 90,923 | 3.64% | 23,027 | 67,896 | 42,542 | 47,174 | 0.510 | 12.83% | nRNA chr2:178,525,988-178,804,642(-) |
| chr1:3 | chr1:35,208-248,936,807 | 39,734 | 1.59% | 14,932 | 24,802 | 24,872 | 11,975 | 0.551 | 11.90% | nRNA chr1:152,168,124-152,332,686(+) |
| chr7:3 | chr7:997,338-158,587,769 | 30,914 | 1.24% | 9,859 | 21,055 | 23,129 | 7,301 | 0.463 | 6.07% | nRNA chr7:141,995,878-142,106,747(+) |
| chr11:3 | chr11:87,377-130,296,836 | 30,852 | 1.23% | 10,125 | 20,727 | 16,551 | 13,969 | 0.544 | 15.10% | nRNA ENST00000698129.1 |
| chr3:3 | chr3:1,092,672-198,225,656 | 26,379 | 1.06% | 9,286 | 17,093 | 18,899 | 7,257 | 0.503 | 7.44% | nRNA chr3:38,039,204-38,124,025(+) |
| chr5:3 | chr5:5,422,781-181,464,195 | 25,682 | 1.03% | 8,170 | 17,512 | 13,351 | 11,421 | 0.585 | 16.39% | nRNA chr5:141,136,682-141,241,744(-) |
| chr9:3 | chr9:2,621,777-138,331,545 | 24,276 | 0.97% | 9,662 | 14,614 | 16,166 | 6,880 | 0.553 | 11.72% | nRNA chr9:110,365,247-110,579,741(-) |
| chr8:3 | chr8:120,387-140,451,378 | 22,956 | 0.92% | 6,885 | 16,071 | 16,637 | 6,136 | 0.503 | 7.79% | nRNA chr8:109,362,460-109,537,207(+) |
| chr12:3 | chr12:44,928-133,204,879 | 19,396 | 0.78% | 4,487 | 14,909 | 15,163 | 3,950 | 0.467 | 4.35% | nRNA chr12:40,393,447-40,570,646(+) |
| chr4:3 | chr4:1,895,797-190,195,896 | 18,805 | 0.75% | 7,570 | 11,235 | 11,949 | 6,650 | 0.513 | 8.88% | nRNA chr4:78,057,322-78,544,269(+) |
| chr17:3 | chr17:137,424-83,209,479 | 18,544 | 0.74% | 3,960 | 14,584 | 9,521 | 8,566 | 0.540 | 7.61% | nRNA chr17:28,614,881-28,645,454(-) |
| chr15:3 | chr15:20,438,654-101,876,462 | 15,122 | 0.60% | 3,885 | 11,237 | 7,719 | 7,274 | 0.520 | 4.81% | nRNA chr15:81,331,216-81,374,213(-) |
| chr6:3 | chr6:95,134-170,745,884 | 15,044 | 0.60% | 5,025 | 10,019 | 11,508 | 3,263 | 0.453 | 5.42% | nRNA chr6:56,457,986-56,642,996(-) |
| chr19:3 | chr19:156,842-58,599,378 | 9,281 | 0.37% | 4,140 | 5,141 | 3,859 | 5,301 | 0.547 | 12.32% | mRNA ZNF615 ENST00000618487.4 |
| chr7:883 | chr7:42,757,557-156,969,973 | 8,382 | 0.34% | 2,516 | 5,866 | 4,027 | 3,873 | 0.510 | 7.60% | nRNA chr7:100,197,473-100,214,387(+) |
| chr5:1380 | chr5:62,306,180-150,396,827 | 7,078 | 0.28% | 916 | 6,162 | 6,670 | 342 | 0.717 | 30.59% | nRNA ENST00000611950.1 |
| chr10:3 | chr10:5,765,702-133,468,290 | 6,340 | 0.25% | 2,404 | 3,936 | 5,718 | 434 | 0.536 | 13.49% | mRNA ZNF518A ENST00000624776.4 |
| chr17:867 | chr17:29,562,314-30,185,883 | 5,529 | 0.22% | 1,540 | 3,989 | 1,471 | 4,042 | 0.525 | 8.72% | nRNA chr17:29,560,546-29,707,090(+) |
| chr16:3 | chr16:1,510,646-90,175,600 | 5,458 | 0.22% | 2,753 | 2,705 | 3,428 | 1,693 | 0.523 | 19.24% | mRNA ZFHX3 ENST00000268489.10 |
| chr13:5941 | chr13:96,453,381-110,512,470 | 5,338 | 0.21% | 283 | 5,055 | 3,475 | 1,782 | 0.438 | 4.05% | nRNA chr13:110,149,060-110,299,118(-) |
| chr11:3225 | chr11:5,248,270-64,321,809 | 4,754 | 0.19% | 1,240 | 3,514 | 2,733 | 2,000 | 0.709 | 36.18% | nRNA chr11:64,318,120-64,321,811(+) |
| chr19:2804 | chr19:36,054,750-55,456,484 | 4,624 | 0.18% | 1,351 | 3,273 | 3,972 | 512 | 0.458 | 5.32% | nRNA chr19:42,325,634-42,378,769(+) |
| chr3:1097 | chr3:49,940,368-190,656,820 | 4,062 | 0.16% | 738 | 3,324 | 3,611 | 427 | 0.486 | 4.73% | nRNA chr3:52,495,337-52,524,495(+) |
| chr14:3 | chr14:16,373,182-83,170,058 | 4,006 | 0.16% | 2,095 | 1,911 | 2,180 | 1,556 | 0.585 | 21.54% | mRNA ADAM20 ENST00000256389.5 |
| chrX:21523 | chrX:67,521,831-67,853,372 | 3,992 | 0.16% | 3,423 | 569 | 3,961 | 1 | 0.371 | 9.32% | mRNA AR ENST00000396043.4 |
| chr11:3926 | chr11:67,351,592-85,771,062 | 3,822 | 0.15% | 1,120 | 2,702 | 2,176 | 1,626 | 0.476 | 5.39% | mRNA SYTL2 ENST00000634661.1 |
| chr18:3 | chr18:1,272,287-65,880,930 | 3,725 | 0.15% | 1,487 | 2,238 | 1,827 | 1,858 | 0.540 | 19.01% | nRNA chr18:58,535,414-58,538,552(+) |
| chrX:3 | chrX:7,049,962-155,568,141 | 3,671 | 0.15% | 1,551 | 2,120 | 2,405 | 1,152 | 0.539 | 8.50% | mRNA FRMPD4 ENST00000657982.1 |
| chr20:3 | chr20:13,849,252-64,292,490 | 3,595 | 0.14% | 1,939 | 1,656 | 2,630 | 907 | 0.615 | 28.82% | mRNA CHD6 ENST00000373233.8 |
| chr17:8543 | chr17:10,300,830-46,371,505 | 3,548 | 0.14% | 123 | 3,425 | 13 | 3,512 | 0.437 | 3.66% | nRNA chr17:10,383,155-10,537,862(+) |

## Top Assigned RNA Targets

| Target | Observed region | False RNA | False frac | mRNA | nRNA | Same | Opp | Mean ZW |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| nRNA chr2:178,525,988-178,804,642(-) | chr2:178,526,810-178,804,640 | 15,139 | 0.61% | 0 | 15,139 | 436 | 14,677 | 0.766 |
| nRNA chr2:178,523,417-178,776,534(+) | chr2:178,526,721-178,776,531 | 7,538 | 0.30% | 0 | 7,538 | 0 | 7,527 | 0.483 |
| nRNA chr11:64,251,529-64,267,923(+) | chr11:64,251,522-64,328,765 | 2,931 | 0.12% | 0 | 2,931 | 2,923 | 0 | 0.767 |
| mRNA ZFHX4 ENST00000518282.5 | chr8:76,704,049-92,818,668 | 1,886 | 0.08% | 1,886 | 0 | 1,885 | 1 | 0.774 |
| nRNA chr5:141,136,682-141,241,744(-) | chr5:141,136,836-141,241,485 | 1,855 | 0.07% | 0 | 1,855 | 2 | 1,832 | 0.683 |
| nRNA chr9:110,365,247-110,579,741(-) | chr9:110,366,265-110,579,724 | 1,851 | 0.07% | 0 | 1,851 | 1,848 | 0 | 0.644 |
| mRNA MUC16 ENST00000397910.8 | chr19:8,886,745-8,981,251 | 1,814 | 0.07% | 1,814 | 0 | 1,807 | 0 | 0.548 |
| nRNA ENST00000698129.1 | chr11:65,499,302-65,506,771 | 1,811 | 0.07% | 0 | 1,811 | 0 | 1,801 | 0.858 |
| mRNA FSIP2 ENST00000424728.6 | chr2:185,788,659-185,833,291 | 1,804 | 0.07% | 1,804 | 0 | 648 | 1,153 | 0.819 |
| nRNA chr1:152,168,124-152,332,686(+) | chr1:152,172,762-152,326,058 | 1,764 | 0.07% | 0 | 1,764 | 13 | 1,292 | 0.795 |
| nRNA chr8:109,362,460-109,537,207(+) | chr8:109,362,514-109,534,214 | 1,757 | 0.07% | 0 | 1,757 | 1,743 | 0 | 0.613 |
| nRNA chr17:29,560,546-29,707,090(+) | chr17:29,562,370-29,703,284 | 1,738 | 0.07% | 0 | 1,738 | 198 | 1,532 | 0.508 |
| nRNA chr21:8,380,664-8,450,668(+) | mixed:101,586-45,024,528 | 1,714 | 0.07% | 0 | 1,714 | 128 | 0 | 0.830 |
| mRNA ZNF462 ENST00000277225.10 | chr9:29,052,401-107,011,284 | 1,606 | 0.06% | 1,606 | 0 | 1,585 | 18 | 0.654 |
| mRNA FLG2 ENST00000388718.5 | chr1:152,350,091-152,357,638 | 1,581 | 0.06% | 1,581 | 0 | 0 | 1,561 | 0.770 |
| mRNA AHNAK ENST00000378024.9 | chr11:62,516,540-62,534,074 | 1,577 | 0.06% | 1,577 | 0 | 1,565 | 0 | 0.732 |
| mRNA AR ENST00000396043.4 | chrX:67,544,227-67,724,281 | 1,563 | 0.06% | 1,563 | 0 | 1,556 | 0 | 0.498 |
| nRNA chr7:141,995,878-142,106,747(+) | chr7:142,005,363-142,106,213 | 1,537 | 0.06% | 0 | 1,537 | 1,513 | 0 | 0.650 |
| nRNA chr4:78,057,322-78,544,269(+) | chr4:78,057,744-78,541,207 | 1,511 | 0.06% | 0 | 1,511 | 1,504 | 0 | 0.501 |
| nRNA chr1:152,335,378-152,364,080(+) | chr1:152,344,441-152,359,090 | 1,505 | 0.06% | 0 | 1,505 | 1 | 1,480 | 0.817 |
| nRNA chr2:151,485,335-151,734,487(-) | chr2:151,485,499-151,733,434 | 1,475 | 0.06% | 0 | 1,475 | 1,332 | 139 | 0.455 |
| nRNA chr11:64,823,051-64,844,653(-) | chr11:64,687,206-64,844,644 | 1,452 | 0.06% | 0 | 1,452 | 1,451 | 0 | 0.683 |
| nRNA chr2:178,528,764-178,620,148(+) | chr2:32,916,265-178,620,133 | 1,450 | 0.06% | 0 | 1,450 | 0 | 1,446 | 0.133 |
| mRNA PCLO ENST00000423517.6 | chr7:82,821,985-83,162,894 | 1,438 | 0.06% | 1,438 | 0 | 1,429 | 0 | 0.687 |
| nRNA ENST00000242310.4 | chr9:32,629,863-32,635,660 | 1,429 | 0.06% | 0 | 1,429 | 0 | 580 | 0.865 |
| mRNA FAT1 ENST00000441802.7 | chr4:83,052,074-186,709,840 | 1,419 | 0.06% | 1,419 | 0 | 1,415 | 1 | 0.807 |
| nRNA chr13:110,149,060-110,299,118(-) | chr13:110,149,064-110,298,934 | 1,415 | 0.06% | 0 | 1,415 | 1,406 | 0 | 0.363 |
| mRNA IRS2 ENST00000375856.5 | chr13:109,755,835-109,786,239 | 1,392 | 0.06% | 1,392 | 0 | 1,215 | 177 | 0.841 |
| nRNA chr6:33,621,321-33,696,574(+) | chr6:33,621,491-33,892,590 | 1,392 | 0.06% | 0 | 1,392 | 1,338 | 45 | 0.530 |
| mRNA SON ENST00000300278.8 | chr21:33,549,478-33,560,569 | 1,362 | 0.05% | 1,362 | 0 | 1,304 | 57 | 0.607 |

## Representative False-Positive Fragments From Top Loci

| Locus | QNAME | Region | Detail | CIGAR | TLEN | NH | NM | ZW | ZC | ZS | Target |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr2:3 | A00839:75:H7MFFDSXY:3:2324:24279:31955 | chr2:202,765,716-202,765,969 | mrna | 150M | -254 | 1 | 0 | 0.9696058630943298 | ambig_same_strand | unspliced | mRNA FAM117B ENST00000392238.3 |
| chr2:3 | A00839:75:H7MFFDSXY:3:1265:27353:27696 | chr2:70,964,870-70,965,128 | nrna | 150M | -259 | 1 | 0 | 0.06471370905637741 | ambig_opp_strand | unspliced | nRNA chr2:70,935,920-70,965,401(+) |
| chr2:3 | A00839:75:H7MFFDSXY:3:2171:22444:6950 | chr2:108,775,607-108,775,927 | mrna | 150M | -321 | 1 | 0 | 0.019881365820765495 | ambig_same_strand | unspliced | mRNA RANBP2 ENST00000697749.1 |
| chr2:3 | A00839:75:H7MFFDSXY:3:2571:30174:27696 | chr2:178,635,925-178,636,234 | nrna | 150M | -310 | 1 | 0 | 0.6275390982627869 | ambig_opp_strand | unspliced | nRNA chr2:178,523,417-178,776,534(+) |
| chr1:3 | A00839:75:H7MFFDSXY:3:1376:18304:15781 | chr1:109,750,676-109,750,856 | nrna | 150M | 181 | 1 | 0 | 0.7252898216247559 | ambig_opp_strand | unspliced | nRNA chr1:109,750,355-109,751,546(-) |
| chr1:3 | A00839:75:H7MFFDSXY:3:1436:25536:12712 | chr1:152,350,515-152,350,871 | mrna | 150M | 357 | 1 | 0 | 0.6577302813529968 | ambig_opp_strand | unspliced | mRNA FLG2 ENST00000388718.5 |
| chr1:3 | A00839:75:H7MFFDSXY:3:2130:1289:21762 | chr1:152,353,387-152,353,563 | nrna | 150M | -177 | 1 | 1 | 0.9367368817329407 | ambig_opp_strand | unspliced | nRNA chr1:152,335,378-152,364,080(+) |
| chr1:3 | A00839:75:H7MFFDSXY:3:2674:3161:10441 | chr1:152,351,717-152,351,998 | nrna | 150M | -282 | 1 | 0 | 0.9113338589668274 | ambig_opp_strand | unspliced | nRNA chr1:152,335,378-152,364,080(+) |
| chr7:3 | A00839:75:H7MFFDSXY:3:2571:8169:34898 | chr7:117,297,519-117,297,759 | nrna | 150M | 241 | 1 | 0 | 0.5559603571891785 | ambig_opp_strand | unspliced | nRNA chr7:117,277,530-117,323,114(-) |
| chr7:3 | A00839:75:H7MFFDSXY:3:1342:15275:5963 | chr7:37,920,560-37,920,725 | nrna | 150M | -166 | 1 | 0 | 0.19380192458629608 | ambig_opp_strand | unspliced | nRNA chr7:37,683,842-37,950,412(+) |
| chr7:3 | A00839:75:H7MFFDSXY:3:1666:13919:21277 | chr7:90,163,010-90,163,350 | nrna | 150M | -341 | 1 | 0 | 0.4697644114494324 | ambig_opp_strand | unspliced | nRNA chr7:90,161,265-90,164,829(+) |
| chr7:3 | A00839:75:H7MFFDSXY:3:1155:15709:16986 | chr7:141,919,494-141,919,677 | nrna | 150M | -184 | 1 | 0 | 0.9467695355415344 | ambig_same_strand | unspliced | nRNA ENST00000548136.1 |
| chr11:3 | A00839:75:H7MFFDSXY:3:2218:30761:29277 | chr11:57,737,560-57,737,844 | nrna | 150M | -285 | 1 | 0 | 0.3247005343437195 | ambig_same_strand | unspliced | nRNA chr11:57,712,602-57,740,369(+) |
| chr11:3 | A00839:75:H7MFFDSXY:3:2224:27778:12054 | chr11:65,223,531-65,223,708 | nrna | 150M | -178 | 1 | 0 | 0.12280038744211197 | ambig_same_strand | unspliced | nRNA chr11:65,213,872-65,242,355(+) |
| chr11:3 | A00839:75:H7MFFDSXY:3:1234:31620:12477 | chr11:64,219,947-64,220,282 | nrna | 150M | -336 | 1 | 0 | 0.39328551292419434 | ambig_same_strand | unspliced | nRNA chr11:64,219,509-64,220,669(+) |
| chr11:3 | A00839:75:H7MFFDSXY:3:2675:25003:35055 | chr11:64,225,500-64,225,802 | nrna | 150M | 303 | 1 | 0 | 0.4868030250072479 | ambig_same_strand | unspliced | nRNA chr11:64,223,798-64,226,158(-) |
| chr3:3 | A00839:75:H7MFFDSXY:3:2437:16893:33786 | chr3:47,343,515-47,343,751 | nrna | 150M | -237 | 1 | 0 | 0.7031286358833313 | ambig_same_strand | unspliced | nRNA chr3:47,282,943-47,344,577(+) |
| chr3:3 | A00839:75:H7MFFDSXY:3:1233:16206:23296 | chr3:53,800,253-53,800,507 | nrna | 150M | -255 | 1 | 0 | 0.7317034006118774 | ambig_same_strand | unspliced | nRNA chr3:53,799,890-53,813,712(+) |
| chr3:3 | A00839:75:H7MFFDSXY:3:2552:23764:25801 | chr3:47,891,475-47,891,772 | mrna | 150M | 298 | 1 | 1 | 0.8284867405891418 | ambig_same_strand | unspliced | mRNA MAP4 ENST00000497735.1 |
| chr3:3 | A00839:75:H7MFFDSXY:3:2132:12274:11412 | chr3:33,514,509-33,514,814 | nrna | 150M | 306 | 1 | 0 | 0.2172052264213562 | ambig_same_strand | unspliced | nRNA chr3:33,498,512-33,645,469(-) |
| chr5:3 | A00839:75:H7MFFDSXY:3:1163:14262:25441 | chr5:141,183,091-141,183,317 | nrna | 150M | 227 | 1 | 0 | 0.870143711566925 | ambig_opp_strand | unspliced | nRNA chr5:141,136,682-141,241,744(-) |
| chr5:3 | A00839:75:H7MFFDSXY:3:1512:7401:2409 | chr5:141,235,450-141,235,769 | nrna | 150M | -320 | 1 | 0 | 0.22802206873893738 | ambig_opp_strand | unspliced | nRNA ENST00000526308.2 |
| chr5:3 | A00839:75:H7MFFDSXY:3:2438:24144:12837 | chr5:140,652,894-140,653,167 | nrna | 150M | -274 | 1 | 0 | 0.6530631184577942 | ambig_same_strand | unspliced | nRNA chr5:140,647,828-140,662,480(+) |
| chr5:3 | A00839:75:H7MFFDSXY:3:1562:25310:29230 | chr5:81,795,907-169,882,931 | mrna | 150M | 88087025 | 1 | 0 | 1.0 | ambig_opp_strand | spliced_implicit | mRNA INSYN2B ENST00000377365.4 |
| chr9:3 | A00839:75:H7MFFDSXY:3:1143:16740:3208 | chr9:120,471,757-120,471,980 | nrna | 150M | 224 | 1 | 0 | 0.36762532591819763 | ambig_same_strand | unspliced | nRNA chr9:120,406,145-120,536,437(-) |
| chr9:3 | A00839:75:H7MFFDSXY:3:2504:2248:3568 | chr9:92,615,947-92,616,184 | mrna | 150M | -238 | 1 | 0 | 0.8769866824150085 | ambig_opp_strand | unspliced | mRNA CENPP ENST00000375587.8 |
| chr9:3 | A00839:75:H7MFFDSXY:3:2512:9552:36010 | chr9:106,924,052-106,924,272 | nrna | 150M | -221 | 1 | 1 | 0.9707114100456238 | ambig_same_strand | unspliced | nRNA chr9:106,863,165-106,930,575(+) |
| chr9:3 | A00839:75:H7MFFDSXY:3:2358:21423:33833 | chr9:110,481,103-110,481,314 | nrna | 150M | 212 | 1 | 0 | 0.8507441282272339 | ambig_same_strand | unspliced | nRNA chr9:110,365,247-110,579,741(-) |
| chr8:3 | A00839:75:H7MFFDSXY:3:1517:4110:30279 | chr8:99,891,808-99,892,048 | nrna | 150M | 241 | 1 | 1 | 0.6409411430358887 | ambig_same_strand | unspliced | nRNA chr8:99,887,060-99,893,697(-) |
| chr8:3 | A00839:75:H7MFFDSXY:3:1343:12310:18959 | chr8:72,937,822-72,937,992 | mrna | 150M | -171 | 1 | 0 | 0.9702807068824768 | ambig_same_strand | unspliced | mRNA KCNB2 ENST00000523207.2 |

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

- Error categories: `results/vcap_gdna_false_rna_hotspots_exon_strand_deconv_v1_2026-05-20/gdna_false_rna_categories.tsv`
- Hotspot windows: `results/vcap_gdna_false_rna_hotspots_exon_strand_deconv_v1_2026-05-20/gdna_false_rna_windows.tsv`
- Hotspot EM loci: `results/vcap_gdna_false_rna_hotspots_exon_strand_deconv_v1_2026-05-20/gdna_false_rna_loci.tsv`
- Assigned RNA targets: `results/vcap_gdna_false_rna_hotspots_exon_strand_deconv_v1_2026-05-20/gdna_false_rna_targets.tsv`
- Representative fragments: `results/vcap_gdna_false_rna_hotspots_exon_strand_deconv_v1_2026-05-20/sample_false_rna_fragments.tsv`
- High-confidence examples: `results/vcap_gdna_false_rna_hotspots_exon_strand_deconv_v1_2026-05-20/high_confidence_false_rna_fragments.tsv`
- Summary JSON: `results/vcap_gdna_false_rna_hotspots_exon_strand_deconv_v1_2026-05-20/summary.json`
