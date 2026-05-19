# VCaP RNA/gDNA Real Mix Pool Confusion - 2026-05-19

BAM: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/with_mm/annotated.bam`

Counting rule: one primary read1 record per fragment; secondary, supplementary, and read2 records are skipped. Truth source is derived from the query-name flowcell ID.

Source mapping:

- `C6EL5ANXX` -> `rna`
- `H7MFFDSXY` -> `gdna`

Fragments counted: 32,218,601

## Headline

Rigel reports 47.88% gDNA and 50.81% RNA across all counted fragments, with 1.31% unresolved. Among assigned RNA/gDNA fragments only, the reported composition is 48.51% gDNA and 51.49% RNA.

Using the flowcell labels as truth, RNA-source fragments are called RNA at 96.21% and gDNA at 2.91%. gDNA-source fragments are called gDNA at 82.40% and RNA at 15.98%.

## Flowcell Summary

| Flowcell | Truth | Fragments | Mean read len | Pred RNA | Pred gDNA | Pred unresolved |
| --- | --- | --- | --- | --- | --- | --- |
| C6EL5ANXX | rna | 13,989,924 | 124.6 | 96.21% | 2.91% | 0.89% |
| H7MFFDSXY | gdna | 18,228,677 | 149.4 | 15.98% | 82.40% | 1.63% |

## RNA versus gDNA Confusion Matrix

| True source | Fragments | Pred RNA | Pred gDNA | Pred unresolved |
| --- | --- | --- | --- | --- |
| rna | 13,989,924 | 13,459,053 (96.21%) | 406,524 (2.91%) | 124,347 (0.89%) |
| gdna | 18,228,677 | 2,912,432 (15.98%) | 15,019,783 (82.40%) | 296,462 (1.63%) |

## Detailed Rigel Pool Calls

| True source | Fragments | Pred mRNA | Pred nRNA | Pred gDNA | Pred unresolved |
| --- | --- | --- | --- | --- | --- |
| rna | 13,989,924 | 12,959,647 (92.64%) | 499,406 (3.57%) | 406,524 (2.91%) | 124,347 (0.89%) |
| gdna | 18,228,677 | 1,340,635 (7.35%) | 1,571,797 (8.62%) | 15,019,783 (82.40%) | 296,462 (1.63%) |

## Largest RNA/gDNA Pool Errors

| Flowcell | Truth | Pred | Detail | ZC | ZS | Count |
| --- | --- | --- | --- | --- | --- | --- |
| H7MFFDSXY | gdna | rna | mrna | ambig_same_strand | unspliced | 1,056,095 |
| H7MFFDSXY | gdna | rna | nrna | ambig_same_strand | unspliced | 1,009,703 |
| H7MFFDSXY | gdna | rna | nrna | ambig_opp_strand | unspliced | 446,849 |
| C6EL5ANXX | rna | gdna | gdna | ambig_same_strand | unspliced | 328,391 |
| H7MFFDSXY | gdna | rna | mrna | ambig_opp_strand | unspliced | 248,193 |
| H7MFFDSXY | gdna | unresolved | unresolved | missing | missing | 147,858 |
| H7MFFDSXY | gdna | unresolved | unresolved | . | unspliced | 145,027 |
| C6EL5ANXX | rna | unresolved | unresolved | missing | missing | 98,133 |
| C6EL5ANXX | rna | gdna | gdna | ambig_opp_strand | unspliced | 54,888 |
| H7MFFDSXY | gdna | rna | nrna | ambig_same_strand | spliced_implicit | 41,963 |
| H7MFFDSXY | gdna | rna | nrna | unambig | unspliced | 32,188 |
| H7MFFDSXY | gdna | rna | nrna | multimapper | unspliced | 15,303 |
| C6EL5ANXX | rna | unresolved | unresolved | . | spliced_annot | 11,966 |
| C6EL5ANXX | rna | unresolved | unresolved | . | unspliced | 11,703 |
| H7MFFDSXY | gdna | rna | nrna | multimapper | spliced_annot | 10,570 |
| C6EL5ANXX | rna | gdna | gdna | unambig | unspliced | 10,374 |
| C6EL5ANXX | rna | gdna | gdna | . | unknown | 9,491 |
| H7MFFDSXY | gdna | rna | mrna | ambig_same_strand | spliced_annot | 9,159 |
| H7MFFDSXY | gdna | rna | mrna | multimapper | unspliced | 7,134 |
| H7MFFDSXY | gdna | rna | nrna | ambig_opp_strand | spliced_implicit | 6,826 |

## Artifacts

- Flowcell summary: `results/vcap_rna20m_gdna20m_with_mm_confusion_2026-05-19/flowcell_summary.tsv`
- RNA/gDNA confusion counts: `results/vcap_rna20m_gdna20m_with_mm_confusion_2026-05-19/confusion_counts.tsv`
- RNA/gDNA confusion matrix: `results/vcap_rna20m_gdna20m_with_mm_confusion_2026-05-19/confusion_matrix.tsv`
- Detailed predicted pools: `results/vcap_rna20m_gdna20m_with_mm_confusion_2026-05-19/detailed_predicted_pool.tsv`
- ZF by flowcell: `results/vcap_rna20m_gdna20m_with_mm_confusion_2026-05-19/zf_by_flowcell.tsv`
- Error ZC/ZS breakdown: `results/vcap_rna20m_gdna20m_with_mm_confusion_2026-05-19/error_zc_zs_breakdown.tsv`
- Summary JSON: `results/vcap_rna20m_gdna20m_with_mm_confusion_2026-05-19/summary.json`
