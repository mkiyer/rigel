# VCaP No-MM gDNA -> RNA Deep Region Analysis - 2026-05-18

BAM: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/no_mm/annotated.bam`

This report focuses on five no-multimap regions with severe gDNA-source fragments miscalled as RNA. Coordinates are 1-based closed for browser use in tables; BED artifacts are 0-based half-open.

Important limitation: the annotated BAM stores the winner (`ZT`/`ZF`) and winner posterior (`ZW`), but not the full candidate log-likelihood vector or runner-up gDNA posterior. The per-fragment likelihood table therefore reconstructs the decision drivers from available run outputs: empirical RNA/gDNA fragment-length likelihood ratio from this BAM, strand likelihood ratio from the run's R1-antisense model, target abundance/effective-length versus locus gDNA abundance/effective-length, and the observed `ZW`. This is enough to identify the mechanism, but adding runner-up gDNA posterior tags would make this exact.

## Selected Regions

| Region | Label | gDNA fragments | gDNA->RNA | False rate | mRNA | nRNA | RNA fragments | Mean ZW | ZW >= .90 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| chr11:65,500,001-65,510,000 | MALAT1 / TALAM1 neighborhood | 4,206 | 3,822 | 90.87% | 1,981 | 1,841 | 64,840 | 0.631 | 36.42% |
| chrX:67,540,001-67,550,000 | AR exon-rich window | 6,585 | 3,270 | 49.66% | 3,256 | 14 | 27,504 | 0.299 | 0.64% |
| chr2:178,530,001-178,620,000 | Long nRNA span around chr2:178.53-178.62 Mb | 26,624 | 18,664 | 70.10% | 1,135 | 17,529 | 1,181 | 0.404 | 9.05% |
| chr1:152,350,001-152,360,000 | FLG2 / FLG repeat-like exon window | 4,610 | 2,868 | 62.21% | 1,503 | 1,365 | 8 | 0.680 | 0.66% |
| chr9:106,920,001-106,930,000 | ZNF462 expressed exon window | 4,798 | 2,839 | 59.17% | 2,689 | 150 | 11,077 | 0.620 | 39.70% |

## Densest 1 kb Sub-Hotspots

These are the first places to open in IGV/browser inside each selected window.

| Region | Densest sub-hotspot | gDNA->RNA fragments | Share of selected-region false calls |
| --- | --- | --- | --- |
| MALAT1 / TALAM1 neighborhood | chr11:65,504,001-65,505,000 | 826 | 21.61% |
| AR exon-rich window | chrX:67,545,001-67,546,000 | 2,172 | 66.42% |
| Long nRNA span around chr2 | chr2:178,539,001-178,540,000 | 363 | 1.94% |
| FLG2 / FLG repeat-like exon window | chr1:152,356,001-152,357,000 | 517 | 18.03% |
| ZNF462 expressed exon window | chr9:106,928,001-106,929,000 | 684 | 24.09% |

## Per-Fragment Examples

The example rows below are chosen from the dominant ambiguous-unspliced false-positive category for each region when available. Columns `fl lr`, `strand lr`, and `abund lr` are log RNA-over-gDNA terms. `p2` is a two-pool RNA probability reconstructed from those terms and should be read as diagnostic, not as the exact Rigel posterior.

## Likelihood-Term Summary

| Region | n | Median ZW | Median fl lr | Median strand lr | Median abund lr | Median p2 | p2 >= .90 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| MALAT1 / TALAM1 neighborhood | 3,758 | 0.700 | -0.290 | 0.692 | 3.728 | 98.64% | 98.88% |
| AR exon-rich window | 3,266 | 0.322 | -0.274 | 0.692 | 4.255 | 98.34% | 76.76% |
| Long nRNA span around chr2:178.53-178.62 Mb | 18,452 | 0.257 | -0.129 | 0.692 | 0.951 | 82.52% | 20.18% |
| FLG2 / FLG repeat-like exon window | 2,847 | 0.792 | -0.098 | 0.692 | 2.118 | 85.37% | 34.98% |
| ZNF462 expressed exon window | 2,839 | 0.553 | -0.205 | 0.692 | 6.129 | 99.69% | 75.80% |

### MALAT1 / TALAM1 neighborhood

Browser region: `chr11:65,500,001-65,510,000`

Top no-MM 10kb window: 3,822 gDNA->RNA calls, 90.9% local false-RNA rate.

| QNAME | Span | Pred | Target | TLEN | CIGARs | NH | NM | ZC | ZS | ZW | fl lr | strand lr | abund lr | p2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A00839:75:H7MFFDSXY:3:2353:32470:32268 | chr11:65,506,277-65,506,501 | nrna | ENST00000698129.1 | 225 | 150M/150M | 1 | 1 | ambig_opp_strand | unspliced | 0.986 | 0.472 | 0.692 | 3.106 | 98.62% |
| A00839:75:H7MFFDSXY:3:1103:31268:6511 | chr11:65,506,478-65,506,701 | nrna | ENST00000698129.1 | 224 | 150M/150M | 1 | 0 | ambig_opp_strand | unspliced | 0.986 | 0.464 | 0.692 | 3.106 | 98.61% |
| A00839:75:H7MFFDSXY:3:1612:2266:27242 | chr11:65,506,380-65,506,603 | nrna | ENST00000698129.1 | 224 | 150M/150M | 1 | 0 | ambig_opp_strand | unspliced | 0.986 | 0.464 | 0.692 | 3.106 | 98.61% |
| A00839:75:H7MFFDSXY:3:1414:8784:21746 | chr11:65,506,459-65,506,671 | nrna | ENST00000698129.1 | 213 | 150M/150M | 1 | 0 | ambig_opp_strand | unspliced | 0.986 | 0.497 | 0.692 | 3.106 | 98.66% |
| A00839:75:H7MFFDSXY:3:1528:5430:14184 | chr11:65,506,279-65,506,491 | nrna | ENST00000698129.1 | 213 | 150M/150M | 1 | 0 | ambig_opp_strand | unspliced | 0.986 | 0.497 | 0.692 | 3.106 | 98.66% |
| A00839:75:H7MFFDSXY:3:2153:5520:35978 | chr11:65,506,345-65,506,579 | nrna | ENST00000698129.1 | 235 | 150M/150M | 1 | 0 | ambig_opp_strand | unspliced | 0.986 | 0.391 | 0.692 | 3.106 | 98.51% |
| A00839:75:H7MFFDSXY:3:2374:11993:6167 | chr11:65,506,489-65,506,735 | nrna | ENST00000698129.1 | 247 | 64S86M/148M | 1 | 6 | ambig_opp_strand | unspliced | 0.985 | 0.317 | 0.692 | 3.106 | 98.39% |
| A00839:75:H7MFFDSXY:3:2273:14362:28745 | chr11:65,506,351-65,506,567 | nrna | ENST00000534336.3 | -217 | 150M/149M1S | 1 | 0 | ambig_opp_strand | unspliced | 0.984 | 0.500 | 0.692 | 2.933 | 98.41% |

### AR exon-rich window

Browser region: `chrX:67,540,001-67,550,000`

Second no-MM 10kb window: 3,270 gDNA->RNA calls, mostly mRNA AR.

| QNAME | Span | Pred | Target | TLEN | CIGARs | NH | NM | ZC | ZS | ZW | fl lr | strand lr | abund lr | p2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A00839:75:H7MFFDSXY:3:1240:8205:35336 | chrX:67,544,562-67,544,806 | mrna | AR | -245 | 150M/150M | 1 | 0 | ambig_same_strand | unspliced | 0.964 | 0.331 | 0.692 | 4.278 | 99.50% |
| A00839:75:H7MFFDSXY:3:1360:3857:14465 | chrX:67,544,234-67,544,428 | mrna | AR | -195 | 150M/129M21S | 1 | 0 | ambig_same_strand | unspliced | 0.964 | 0.459 | 0.692 | 4.278 | 99.56% |
| A00839:75:H7MFFDSXY:3:1464:21811:23703 | chrX:67,544,503-67,544,751 | mrna | AR | -249 | 150M/150M | 1 | 0 | ambig_same_strand | unspliced | 0.964 | 0.304 | 0.692 | 4.278 | 99.49% |
| A00839:75:H7MFFDSXY:3:2326:7374:13416 | chrX:67,544,708-67,544,962 | mrna | AR | -255 | 150M/150M | 1 | 0 | ambig_same_strand | unspliced | 0.964 | 0.231 | 0.692 | 4.278 | 99.45% |
| A00839:75:H7MFFDSXY:3:2568:19253:22153 | chrX:67,544,790-67,545,063 | mrna | AR | -274 | 150M/150M | 1 | 1 | ambig_same_strand | unspliced | 0.963 | -0.028 | 0.692 | 4.278 | 99.29% |
| A00839:75:H7MFFDSXY:3:2474:31349:17832 | chrX:67,544,257-67,544,534 | mrna | AR | -278 | 150M/93M57S | 1 | 0 | ambig_same_strand | unspliced | 0.963 | -0.088 | 0.692 | 4.278 | 99.25% |
| A00839:75:H7MFFDSXY:3:1373:11469:14497 | chrX:67,544,747-67,545,037 | mrna | AR | -291 | 150M/150M | 1 | 0 | ambig_same_strand | unspliced | 0.961 | -0.270 | 0.692 | 4.278 | 99.10% |
| A00839:75:H7MFFDSXY:3:1102:10176:32111 | chrX:67,544,252-67,544,544 | mrna | AR | -293 | 150M/111M39S | 1 | 0 | ambig_same_strand | unspliced | 0.961 | -0.297 | 0.692 | 4.278 | 99.07% |

### Long nRNA span around chr2:178.53-178.62 Mb

Browser region: `chr2:178,530,001-178,620,000`

Contiguous severe nRNA false-positive band with little local RNA support.

| QNAME | Span | Pred | Target | TLEN | CIGARs | NH | NM | ZC | ZS | ZW | fl lr | strand lr | abund lr | p2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A00839:75:H7MFFDSXY:3:2655:10529:23578 | chr2:178,606,895-178,607,117 | nrna | RIGEL_NRNA_chr2_2_178525988_178804642 | 223 | 150M/150M | 1 | 0 | ambig_opp_strand | unspliced | 0.935 | 0.476 | 0.692 | 1.485 | 93.42% |
| A00839:75:H7MFFDSXY:3:1341:8160:13777 | chr2:178,601,362-178,601,578 | nrna | RIGEL_NRNA_chr2_2_178525988_178804642 | 217 | 150M/150M | 1 | 0 | ambig_opp_strand | unspliced | 0.935 | 0.500 | 0.692 | 1.485 | 93.57% |
| A00839:75:H7MFFDSXY:3:2368:1588:18114 | chr2:178,608,083-178,608,301 | nrna | RIGEL_NRNA_chr2_2_178525988_178804642 | 219 | 150M/150M | 1 | 1 | ambig_opp_strand | unspliced | 0.935 | 0.485 | 0.692 | 1.485 | 93.47% |
| A00839:75:H7MFFDSXY:3:1166:32362:9721 | chr2:178,607,698-178,607,916 | nrna | RIGEL_NRNA_chr2_2_178525988_178804642 | 219 | 150M/150M | 1 | 0 | ambig_opp_strand | unspliced | 0.935 | 0.485 | 0.692 | 1.485 | 93.47% |
| A00839:75:H7MFFDSXY:3:1327:20157:32675 | chr2:178,605,196-178,605,410 | nrna | RIGEL_NRNA_chr2_2_178525988_178804642 | 215 | 150M/150M | 1 | 0 | ambig_opp_strand | unspliced | 0.935 | 0.491 | 0.692 | 1.485 | 93.51% |
| A00839:75:H7MFFDSXY:3:2104:6985:26083 | chr2:178,607,310-178,607,524 | nrna | RIGEL_NRNA_chr2_2_178525988_178804642 | 215 | 150M/150M | 1 | 0 | ambig_opp_strand | unspliced | 0.935 | 0.491 | 0.692 | 1.485 | 93.51% |
| A00839:75:H7MFFDSXY:3:1206:7997:1407 | chr2:178,602,118-178,602,332 | nrna | RIGEL_NRNA_chr2_2_178525988_178804642 | 215 | 150M/150M | 1 | 0 | ambig_opp_strand | unspliced | 0.935 | 0.491 | 0.692 | 1.485 | 93.51% |
| A00839:75:H7MFFDSXY:3:1648:23113:27555 | chr2:178,597,889-178,598,103 | nrna | RIGEL_NRNA_chr2_2_178525988_178804642 | 215 | 150M/150M | 1 | 0 | ambig_opp_strand | unspliced | 0.935 | 0.491 | 0.692 | 1.485 | 93.51% |

### FLG2 / FLG repeat-like exon window

Browser region: `chr1:152,350,001-152,360,000`

High false-RNA rate with nearly no RNA-source fragments in the window.

| QNAME | Span | Pred | Target | TLEN | CIGARs | NH | NM | ZC | ZS | ZW | fl lr | strand lr | abund lr | p2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A00839:75:H7MFFDSXY:3:1255:10465:15608 | chr1:152,352,772-152,352,828 | mrna | FLG2 | 57 | 57M/57M | 1 | 0 | ambig_opp_strand | unspliced | 0.893 | 2.122 | 0.692 | 2.118 | 99.28% |
| A00839:75:H7MFFDSXY:3:2251:5168:34491 | chr1:152,357,185-152,357,407 | mrna | FLG2 | 223 | 150M/150M | 1 | 0 | ambig_opp_strand | unspliced | 0.877 | 0.476 | 0.692 | 2.118 | 96.40% |
| A00839:75:H7MFFDSXY:3:1443:25418:32362 | chr1:152,356,449-152,356,671 | mrna | FLG2 | 223 | 150M/150M | 1 | 0 | ambig_opp_strand | unspliced | 0.877 | 0.476 | 0.692 | 2.118 | 96.40% |
| A00839:75:H7MFFDSXY:3:1166:26223:33849 | chr1:152,356,627-152,356,849 | mrna | FLG2 | 223 | 150M/150M | 1 | 0 | ambig_opp_strand | unspliced | 0.877 | 0.476 | 0.692 | 2.118 | 96.40% |
| A00839:75:H7MFFDSXY:3:2419:29125:7529 | chr1:152,353,255-152,353,477 | mrna | FLG2 | 223 | 150M/150M | 1 | 0 | ambig_opp_strand | unspliced | 0.877 | 0.476 | 0.692 | 2.118 | 96.40% |
| A00839:75:H7MFFDSXY:3:1221:22905:1110 | chr1:152,357,326-152,357,548 | mrna | FLG2 | 223 | 150M/150M | 1 | 0 | ambig_opp_strand | unspliced | 0.877 | 0.476 | 0.692 | 2.118 | 96.40% |
| A00839:75:H7MFFDSXY:3:1209:32145:34209 | chr1:152,354,270-152,354,492 | mrna | FLG2 | 223 | 150M/145M | 1 | 0 | ambig_opp_strand | unspliced | 0.877 | 0.476 | 0.692 | 2.118 | 96.40% |
| A00839:75:H7MFFDSXY:3:1356:21585:36933 | chr1:152,355,483-152,355,699 | mrna | FLG2 | 217 | 150M/150M | 1 | 1 | ambig_opp_strand | unspliced | 0.877 | 0.500 | 0.692 | 2.118 | 96.48% |

### ZNF462 expressed exon window

Browser region: `chr9:106,920,001-106,930,000`

High mRNA false-positive count in a locally expressed gene.

| QNAME | Span | Pred | Target | TLEN | CIGARs | NH | NM | ZC | ZS | ZW | fl lr | strand lr | abund lr | p2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A00839:75:H7MFFDSXY:3:1168:20057:21042 | chr9:106,925,271-106,925,493 | mrna | ZNF462 | -223 | 150M/150M | 1 | 0 | ambig_same_strand | unspliced | 0.998 | 0.476 | 0.692 | 6.129 | 99.93% |
| A00839:75:H7MFFDSXY:3:1628:1090:7075 | chr9:106,926,383-106,926,605 | mrna | ZNF462 | -223 | 150M/150M | 1 | 1 | ambig_same_strand | unspliced | 0.998 | 0.476 | 0.692 | 6.129 | 99.93% |
| A00839:75:H7MFFDSXY:3:1566:5918:35289 | chr9:106,925,204-106,925,426 | mrna | ZNF462 | -223 | 150M/150M | 1 | 0 | ambig_same_strand | unspliced | 0.998 | 0.476 | 0.692 | 6.129 | 99.93% |
| A00839:75:H7MFFDSXY:3:2256:9462:25833 | chr9:106,926,686-106,926,908 | mrna | ZNF462 | -223 | 150M/150M | 1 | 0 | ambig_same_strand | unspliced | 0.998 | 0.476 | 0.692 | 6.129 | 99.93% |
| A00839:75:H7MFFDSXY:3:2140:24089:14403 | chr9:106,924,971-106,925,187 | mrna | ZNF462 | -217 | 150M/150M | 1 | 2 | ambig_same_strand | unspliced | 0.998 | 0.500 | 0.692 | 6.129 | 99.93% |
| A00839:75:H7MFFDSXY:3:1525:13286:33927 | chr9:106,926,144-106,926,360 | mrna | ZNF462 | -217 | 150M/150M | 1 | 0 | ambig_same_strand | unspliced | 0.998 | 0.500 | 0.692 | 6.129 | 99.93% |
| A00839:75:H7MFFDSXY:3:2445:28890:36902 | chr9:106,925,787-106,926,003 | mrna | ZNF462 | -217 | 150M/150M | 1 | 0 | ambig_same_strand | unspliced | 0.998 | 0.500 | 0.692 | 6.129 | 99.93% |
| A00839:75:H7MFFDSXY:3:1636:17517:30984 | chr9:106,924,911-106,925,129 | mrna | ZNF462 | -219 | 150M/150M | 1 | 1 | ambig_same_strand | unspliced | 0.998 | 0.485 | 0.692 | 6.129 | 99.93% |

## Mechanism

1. The dominant false-positive class is not multimapping. The examples above are `NH=1`, usually `150M/150M`, low `NM`, and `ZS=unspliced`. They are perfectly plausible gDNA fragments that also overlap an RNA candidate.
2. The RNA candidate gets a small but systematic strand advantage: about +0.69 log units over gDNA because the library is highly stranded while the gDNA likelihood is strand-symmetric (`log(0.5)`). This is appropriate for true RNA, but it becomes a bias toward RNA when the fragment has no splice or RNA-unique evidence.
3. Fragment length is not the main separator in these cases. The median empirical RNA-over-gDNA fragment-length term is near zero or slightly negative in every selected region, so most of these fragments are not being rescued by an RNA-specific length signature.
4. The decisive term is abundance per effective length. AR and ZNF462 have median abundance/effective-length advantages of +4.25 and +6.13 log units over locus gDNA; MALAT1/TALAM1 is +3.73. Once an ambiguous unspliced fragment enters the EM equivalence class, the local RNA component mass overwhelms the gDNA component even though the fragment itself is non-diagnostic.
5. The low-expression regions show the more concerning form: FLG2 has only 8 RNA-source fragments in the 10 kb window but still calls 2,868 gDNA-source fragments as RNA, and the chr2 long-nRNA region has 18,664 false RNA calls with only 1,181 RNA-source fragments across 90 kb. This looks like long-span nRNA/self-seeding behavior: gDNA-compatible unspliced reads create enough nRNA mass for the EM to keep pulling more such reads into RNA.
6. `ZW` is winner-component posterior, not pool-level RNA posterior. AR has low median `ZW` because the RNA pool is split across AR isoforms, while the reconstructed two-pool RNA-vs-gDNA probability is still high. For debugging and diagnostics we need runner-up gDNA posterior or pool-level posterior tags, not only the winning component posterior.

## Fix Direction

1. For diagnostic calling, treat `ZS=unspliced` plus `ZC=ambig_same_strand/ambig_opp_strand` as weak RNA evidence. Require a positive splice/unique-exon signal or a large RNA-over-gDNA margin before a hard RNA call.
2. Add a conservative hard-call policy that defaults ambiguous unspliced RNA/gDNA competition to gDNA unless runner-up gDNA posterior is safely low.
3. Emit runner-up pool posterior or `log_posterior_rna - log_posterior_gdna` in annotated BAMs. `ZW` tells us the winner can be confident, but not how close gDNA was.
4. Add local guards for long nRNA spans: if a region has little local RNA-source evidence and many gDNA-compatible unspliced reads, prevent nRNA self-seeding from becoming strong positive RNA evidence.

## Artifacts

- Region summary: `results/vcap_no_mm_deep_regions_2026-05-18/selected_regions_summary.tsv`
- Region BED: `results/vcap_no_mm_deep_regions_2026-05-18/selected_regions.bed`
- 1 kb sub-hotspots: `results/vcap_no_mm_deep_regions_2026-05-18/subhotspots_1kb.tsv`
- Per-fragment examples: `results/vcap_no_mm_deep_regions_2026-05-18/per_fragment_examples.tsv`
- Likelihood-term summary: `results/vcap_no_mm_deep_regions_2026-05-18/likelihood_term_summary.tsv`
- Category counts: `results/vcap_no_mm_deep_regions_2026-05-18/region_false_categories.tsv`
- Target counts: `results/vcap_no_mm_deep_regions_2026-05-18/region_false_targets.tsv`
- Locus counts: `results/vcap_no_mm_deep_regions_2026-05-18/region_false_loci.tsv`
