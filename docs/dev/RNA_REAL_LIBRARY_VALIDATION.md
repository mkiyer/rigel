# Real-library validation: frozen count candidate

*2026-10-07. Four serial paired whole-genome runs. No density-reader integration or production source changes.*

The candidate is the same factory-free count prototype used in the panel report: observed intron strand claims and exact conditional strand messages. The existing capture reader remains active. All eight pipelines finished with finite, nonnegative transcript counts and region DNA counts. Paired scan totals, configuration, input identity and quantification hashes were checked.

## DNA fraction

These are inferred DNA fractions (%), not transcript error. Only VCaP has an available mixture truth: **25.18%** by read name.

| Library | 0.7.1 archived | Current | Candidate |
|---|---:|---:|---:|
| lbx0190 | 13.04 | 8.34 | 8.50 |
| lbx0588 | 96.62 | 91.14 | 92.65 |
| mo3021 | 22.48 | 15.46 | 15.87 |
| vcap | 26.31 | 23.44 | 23.95 |

VCaP moves toward the known mixture, but still under-calls DNA. The archived release used its compatible index and scanner; its assignment denominator differs from the current pair. It is a historical reference, not a perfectly matched mechanism contrast. No unlabelled library supplies a transcript-accuracy verdict.

## RNA output changes

Absolute changes sum the per-transcript or per-gene count differences and divide by the current transcript total. These measure output stability, **not accuracy**. Signed total change can be smaller because increases and decreases cancel.

| Library | Transcript absolute change (%) | Gene absolute change (%) | Transcript total change (%) |
|---|---:|---:|---:|
| lbx0190 | 0.133 | 0.072 | -0.066 |
| lbx0588 | 14.179 | 6.030 | -2.614 |
| mo3021 | 0.144 | 0.063 | -0.047 |
| vcap | 4.033 | 0.280 | +0.133 |

LBX0588 redistributes 14.18% of transcript counts and 6.03% of gene counts; half its transcript change lies in ten genes. VCaP redistributes 4.03% of transcript counts but only 0.28% of gene counts. These changes must not be described as uniformly small. The largest DNA shift is in LBX0588, with most of the displaced mass coming from the unspliced-RNA pool. VCaP transcript total changes by approximately 0.13%; the principal pool change is also less unspliced RNA. This is consistent with removing the intron constraint, but does not prove an improvement in the absence of origin truth for those assignments.

## Reproduction and limitations

Receipts are in `.cache/rigel_runs/2026-10-07_release_followup/real/`: eight JSON results and quantification tables, pipeline logs, `manifest.json` and `comparison.json`. `report_real.py` verifies inputs/settings and renders this document. The first VCaP candidate attempt was interrupted during scanning; its process was confirmed absent before a single retry. The incomplete log is preserved and the retry is recorded.

Each pipeline used one scan thread, fractional assignment, seed zero, default calibration/EM threads and `OMP_NUM_THREADS=1`. Exactly one whole-genome pipeline ran at a time. No whole-genome calibration debug dump was made. Load changed materially during VCaP; elapsed times are receipts, not a speed comparison.

The count candidate now has panel and real-library validation. The subsequent
[readiness audit](RNA_COUNT_READINESS_CHECKPOINT.md) completed the suite, preflight and
golden review. On 2026-10-09 the count foundation was integrated and its eight accepted
goldens updated; the complete cleanup chain and integrated package preserve the LBX0190
pipeline bit-for-bit against the original validated candidate. All four available real
libraries are strongly stranded, so genuine unstranded real-data validation remains a gap.
The separate detector-free reader remains outstanding. See
[RNA_FULL_PANEL_VALIDATION.md](RNA_FULL_PANEL_VALIDATION.md) for every simulation stratum.
