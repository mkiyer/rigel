# VCaP EXON Strand-Deconvolution Validation - 2026-05-20

Input BAM: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam`

Fresh run: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/exon_strand_deconv_v1`

Truth source is flowcell-derived: RNA=`C6EL5ANXX`, gDNA=`H7MFFDSXY`. Counts are one primary read1 per fragment from Rigel annotated BAMs.

## Headline

Relative to `kappa_units_fix`, EXON strand deconvolution changed true-gDNA recall from 85.24% to 84.66%. True gDNA called RNA changed by +106,732 fragments; true RNA called gDNA changed by -41,781 fragments.

## Confusion Summary

| Run | gDNA recall | gDNA -> RNA | gDNA -> mRNA | gDNA -> nRNA | RNA -> gDNA | EM gDNA | EM mRNA | EM nRNA |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| v4_3_with_mm | 77.94% | 3,725,095 (20.44%) | 1,346,710 (7.39%) | 2,378,385 (13.05%) | 383,809 (2.74%) | 14,299,329 | 14,716,943 | 2,499,655 |
| kappa_units_fix | 85.24% | 2,393,558 (13.13%) | 1,062,195 (5.83%) | 1,331,363 (7.30%) | 508,599 (3.64%) | 15,755,863 | 14,330,067 | 1,429,997 |
| exon_strand_deconv_v1 | 84.66% | 2,500,290 (13.72%) | 1,171,473 (6.43%) | 1,328,817 (7.29%) | 466,818 (3.34%) | 15,607,342 | 14,505,283 | 1,403,302 |

## Deltas

| Comparison | Δ gDNA recall pp | Δ gDNA -> RNA | Δ gDNA -> mRNA | Δ gDNA -> nRNA | Δ RNA -> gDNA | Δ EM gDNA |
| --- | --- | --- | --- | --- | --- | --- |
| v4_3_with_mm - kappa_units_fix | -7.30 | +1,331,537 | +284,515 | +1,047,022 | -124,790 | -1,456,534 |
| exon_strand_deconv_v1 - kappa_units_fix | -0.59 | +106,732 | +109,278 | -2,546 | -41,781 | -148,521 |

## Calibration And Exposure

| Run | Mode | rho_ref | spread | null | floor rows | neg floors | median weight | rho EI | rho EC | rho composite |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| v4_3_with_mm | regional | 0.0002132 | 5.819 | 0.5358 | 4,364 | 0 | 0.08417 | 0.0299 | 0 | 0 |
| kappa_units_fix | regional | 0.0008764 | 3.484 | 0.3407 | 6 | 0 | 0.6177 | 0.0299 | 0 | 0 |
| exon_strand_deconv_v1 | regional | 0.0003609 | 2.747 | 0.4048 | 1 | 0 | 0.7849 | 0.0299 | 0.01977 | 0.02381 |

## Interpretation

Relative to the older `v4_3_with_mm` run, the current code is much better: true gDNA called RNA is down by 1,224,805 fragments, with gDNA recall up from 77.94% to 84.66%. That gain is still mostly the earlier kappa-units fix, not the new EXON-contained deconvolution itself.

Relative to the immediate `kappa_units_fix` baseline, the result is mixed and not an overall improvement for this VCaP mixture. The new feature reduces true RNA called gDNA by 41,781 fragments, but increases true gDNA called RNA by 106,732 fragments. The regression is almost entirely mRNA: gDNA -> mRNA rises by 109,278, while gDNA -> nRNA changes by only -2,546.

The feature is doing the intended mathematical operation: EXON-contained strand deconvolution estimates a lower exon gDNA density than the boundary channel (`rho_ec=0.01977` versus `rho_ei=0.02990`) and fuses to a lower composite (`rho_comp=0.02381`). That lowers `mean_pi_gdna` from 0.3163 to 0.2191 and lowers summed locus gDNA priors by 576,657 prior-count units. The EM then moves 148,521 fragments out of gDNA overall, mostly into mRNA.

## Root Cause

The remaining and newly increased errors are not primarily an nRNA siphon. Category deltas versus `kappa_units_fix` show the gDNA false-RNA increase is concentrated in unspliced ambiguous exon-compatible reads:

| Error class | kappa_units_fix | exon_strand_deconv_v1 | Delta |
| --- | ---: | ---: | ---: |
| true gDNA -> mRNA, ambig_same_strand, unspliced | 830,287 | 918,169 | +87,882 |
| true gDNA -> mRNA, ambig_opp_strand, unspliced | 197,129 | 217,716 | +20,587 |
| true gDNA -> nRNA, ambig_opp_strand, unspliced | 342,908 | 361,080 | +18,172 |
| true gDNA -> nRNA, ambig_same_strand, unspliced | 877,209 | 854,615 | -22,594 |

This is the expected failure mode when an exon-contained strand correction is applied to a hybrid-capture/high-expression mixture. The contained EXON channel sees a large strand-asymmetric RNA signal and subtracts it, producing a lower gDNA density. That protects some true RNA from being overcalled as gDNA, but it also weakens the prior for true captured gDNA fragments that land inside mature exons where mRNA and gDNA evidence are locally indistinguishable.

The locus-level deltas are diffuse rather than one catastrophic locus. Compared with `kappa_units_fix`, 9,807 loci lose gDNA EM mass, 5,786 gain, and 8,459 are unchanged. The largest single loss is locus `3`, with -15,209 gDNA and +23,512 mRNA, but the total gDNA EM loss is -148,521, so most of the shift is spread across many exonic/high-expression loci.

The hotspot scan confirms the same anatomy. In the fresh run, 96.12% of gDNA -> RNA false positives are `ZS=unspliced`, and 96.97% are either `ambig_same_strand` or `ambig_opp_strand`. The top 10 EM loci account for 13.20% of all false RNA calls, while the top 10 genomic 10 kb windows account for only 1.13%, so the problem is recurrent across many ambiguous genic regions rather than a single coordinate artifact.

Top windows remain biologically interpretable capture/expression hotspots. Examples include `chr11:65,500,001-65,510,000` with 3,809 / 4,206 true-gDNA fragments called RNA, `chr1:152,350,001-152,360,000` with 3,297 / 4,610, and the AR window `chrX:67,540,001-67,550,000` with 2,890 / 6,585. The AR window is still mostly mRNA false positives: 2,797 mRNA and 93 nRNA.

## Recommendations

Do not merge this feature as a default-on replacement for the kappa-units baseline based on this VCaP result. It is useful evidence, but in this assay it makes the global exon prior too weak for true gDNA in captured exons.

The next implementation step should keep EXON-contained deconvolution as a diagnostic or bounded component rather than allowing it to lower exon gDNA density freely. Concrete options worth testing:

1. Use `max(EXON-INTRON, EXON-COMPOSITE)` for the mature-exon gDNA prior in capture-like libraries, while still reporting EXON-contained diagnostics.
2. Add a capture/local-context floor so exon-contained density cannot fall below a robust regional boundary density in high-exposure exon neighborhoods.
3. Split mature-exon gDNA prior construction into boundary-supported and contained-only terms, with contained-only shrinkage toward boundary density when local RNA abundance is high.
4. Add a post-EM or prior diagnostic that flags loci where lowering exon gDNA prior increases mRNA assignment for true-gDNA-like unspliced ambiguous fragments.

The most promising direction is option 2 or 3: use EXON-contained strand deconvolution to detect RNA contamination, but constrain it with local capture exposure so highly expressed/captured mature exons do not become under-protected against gDNA.

## Largest Locus-Level gDNA Gains Versus kappa_units_fix

| locus_id | n EM | Δ gDNA | Δ mRNA | Δ nRNA | prior old | prior new | w old | w new |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 21523 | 59,500 | +195 | +377 | -572 | 536.9 | 474.8 | 0.2304 | 0.3026 |
| 6019 | 13,544 | +164 | +73 | -237 | 973.5 | 919 | 0.2975 | 0.4316 |
| 5 | 10,753 | +144 | +66 | -210 | 1.599e+04 | 1.55e+04 | 0.5374 | 0.6344 |
| 6945 | 7,436 | +131 | +33 | -164 | 1.249e+04 | 1.179e+04 | 0.3135 | 0.3983 |
| 11351 | 7,153 | +125 | +19 | -144 | 2113 | 2057 | 0.368 | 0.5132 |
| 16275 | 4,030 | +106 | +10 | -116 | 1176 | 1129 | 0.6125 | 0.7494 |
| 13259 | 14,828 | +101 | +3 | -104 | 963.9 | 920.4 | 0.3671 | 0.5063 |
| 20671 | 5,677 | +99 | +17 | -116 | 2098 | 1965 | 0.2347 | 0.3493 |
| 17222 | 4,955 | +98 | +19 | -117 | 1963 | 1931 | 0.4061 | 0.5601 |
| 12551 | 3,652 | +95 | +2 | -97 | 509.5 | 455.5 | 0.5607 | 0.7277 |
| 9609 | 11,402 | +94 | -8 | -86 | 1585 | 1568 | 0.3158 | 0.4853 |
| 4518 | 5,264 | +92 | +20 | -112 | 2489 | 2397 | 0.3822 | 0.5351 |

## Largest Locus-Level gDNA Losses Versus kappa_units_fix

| locus_id | n EM | Δ gDNA | Δ mRNA | Δ nRNA | prior old | prior new | w old | w new |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3 | 2,823,632 | -15,209 | +23,512 | -8,303 | 1.06e+06 | 9.934e+05 | 0.1752 | 0.2612 |
| 5941 | 31,819 | -945 | +74 | +871 | 3199 | 3029 | 0.1618 | 0.2249 |
| 1380 | 44,850 | -764 | +711 | +53 | 9590 | 8754 | 0.2927 | 0.402 |
| 247 | 57,164 | -638 | +281 | +357 | 1.276e+04 | 1.196e+04 | 0.1968 | 0.274 |
| 5019 | 17,656 | -553 | +31 | +522 | 3212 | 3138 | 0.3047 | 0.4119 |
| 3914 | 12,969 | -493 | +485 | +8 | 1443 | 1333 | 0.1033 | 0.1586 |
| 1097 | 41,357 | -397 | +103 | +294 | 7366 | 7006 | 0.3594 | 0.4672 |
| 122 | 17,209 | -384 | +119 | +265 | 4651 | 4446 | 0.2465 | 0.3625 |
| 2561 | 23,820 | -379 | +251 | +128 | 6909 | 6484 | 0.1071 | 0.1653 |
| 16782 | 12,024 | -377 | +378 | -1 | 2447 | 2036 | 0.4184 | 0.4978 |
| 3926 | 21,743 | -362 | +27 | +335 | 1879 | 1745 | 0.4879 | 0.5963 |
| 18563 | 11,192 | -362 | +324 | +38 | 1694 | 1567 | 0.1459 | 0.2171 |

## Artifacts

- Run metrics: `results/vcap_exon_strand_deconv_validation_2026-05-20/run_metrics.tsv`
- Detailed confusion comparison: `results/vcap_exon_strand_deconv_validation_2026-05-20/confusion_detail_compare.tsv`
- Locus deltas: `results/vcap_exon_strand_deconv_validation_2026-05-20/all_locus_delta.tsv`
- Top gDNA gains/losses: `results/vcap_exon_strand_deconv_validation_2026-05-20/top_gdna_gain.tsv`, `results/vcap_exon_strand_deconv_validation_2026-05-20/top_gdna_loss.tsv`
- Fresh confusion report: `docs/benchmarks/vcap_rna20m_gdna20m_exon_strand_deconv_v1_confusion_2026-05-20.md`
- Fresh hotspot report: `docs/benchmarks/vcap_gdna_false_rna_hotspots_exon_strand_deconv_v1_2026-05-20.md`
