# VCaP Global High-Tail Exposure Validation - 2026-05-21

Input run family: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m`

Compared runs:

- Q95 baseline: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/exon_strand_deconv_v1`
- Q99.9: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/global_q999_v1`
- Q99.95: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/global_q9995_v1`

Truth source is flowcell-derived: RNA=`C6EL5ANXX`, gDNA=`H7MFFDSXY`. Counts are one primary read1 per fragment from Rigel annotated BAMs.

## Headline

The best high-tail run, `Q99.9`, increases gDNA-to-RNA leakage: true-gDNA fragments called RNA changed from 2,500,290 to 2,894,671, a 15.77% increase relative to Q95.

## Confusion Summary

| Run | rho_ref | gDNA recall | gDNA -> RNA | gDNA -> mRNA | gDNA -> nRNA | RNA -> gDNA | EM gDNA | EM mRNA | EM nRNA |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Q95 | 0.0003609 | 84.66% | 2,500,290 (13.72%) | 1,171,473 (6.43%) | 1,328,817 (7.29%) | 466,818 (3.34%) | 15,607,342 | 14,505,283 | 1,403,302 |
| Q99.9 | 0.006698 | 82.49% | 2,894,671 (15.88%) | 1,085,039 (5.95%) | 1,809,632 (9.93%) | 459,112 (3.28%) | 15,205,189 | 14,371,304 | 1,939,434 |
| Q99.95 | 0.0101 | 82.31% | 2,928,428 (16.06%) | 1,087,668 (5.97%) | 1,840,760 (10.10%) | 455,500 (3.26%) | 15,167,797 | 14,373,188 | 1,974,942 |

## Deltas Versus Q95 Baseline

| Comparison | Delta gDNA recall pp | Delta gDNA -> RNA | Leakage reduction | Delta gDNA -> mRNA | Delta gDNA -> nRNA | Delta RNA -> gDNA | Delta EM gDNA |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Q99.9 - Q95 | -2.16 | +394,381 | -15.77% | -86,434 | +480,815 | -7,706 | -402,153 |
| Q99.95 - Q95 | -2.35 | +428,138 | -17.12% | -83,805 | +511,943 | -11,318 | -439,545 |

## Locus 3 / FLG2 Mega-Locus

| Run | Exposure weight | L_gDNA | Prior EM | EM gDNA | EM mRNA | EM nRNA |
| --- | --- | --- | --- | --- | --- | --- |
| Q95 | 0.261212 | 92,407,095 | 259,499 | 1,446,008 | 1,162,014 | 288,677 |
| Q99.9 | 0.0268825 | 9,510,006 | 26,706 | 1,378,282 | 1,140,289 | 378,128 |
| Q99.95 | 0.0184693 | 6,533,723 | 18,348 | 1,372,478 | 1,139,974 | 384,247 |

## Interpretation

- Q99.9 increased total true-gDNA -> RNA calls by +394,381 fragments versus Q95. The split was -86,434 mRNA and +480,815 nRNA, with RNA -> gDNA changing by -7,706.
- Q99.95 increased total true-gDNA -> RNA calls by +428,138 fragments versus Q95. The split was -83,805 mRNA and +511,943 nRNA, with RNA -> gDNA changing by -11,318.

The high-tail global reference quantile is therefore not sufficient as a standalone fix for the dominant VCaP leakage mode. The key diagnostic is the mRNA/nRNA split: if mRNA leakage falls while nRNA leakage rises, the higher reference scale is mainly changing the RNA sub-state balance rather than removing the false-RNA assignment pressure.

## Artifacts

- Metrics table: `results/vcap_global_high_tail_exposure_2026-05-21/run_metrics.tsv`
- Confusion detail comparison: `results/vcap_global_high_tail_exposure_2026-05-21/confusion_detail_compare.tsv`
