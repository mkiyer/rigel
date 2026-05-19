# VCaP RNA/gDNA Real Mix Pool Confusion - 2026-05-15

BAM: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam`

Counting rule: one primary read1 record per fragment; secondary, supplementary, and read2 records are skipped. Truth source is derived from the query-name flowcell ID.

Source mapping:

- `C6EL5ANXX` -> `rna`
- `H7MFFDSXY` -> `gdna`

Fragments counted: 32,218,601

## Headline

Rigel reports 46.76% gDNA and 52.19% RNA across all counted fragments, with 1.06% unresolved. Among assigned RNA/gDNA fragments only, the reported composition is 47.25% gDNA and 52.75% RNA.

Using the flowcell labels as truth, RNA-source fragments are called RNA at 96.38% and gDNA at 2.76%. gDNA-source fragments are called gDNA at 80.52% and RNA at 18.28%.

## Flowcell Summary

| Flowcell | Truth | Fragments | Mean read len | Pred RNA | Pred gDNA | Pred unresolved |
| --- | --- | --- | --- | --- | --- | --- |
| C6EL5ANXX | rna | 13,989,924 | 124.6 | 96.38% | 2.76% | 0.86% |
| H7MFFDSXY | gdna | 18,228,677 | 149.4 | 18.28% | 80.52% | 1.20% |

## RNA versus gDNA Confusion Matrix

| True source | Fragments | Pred RNA | Pred gDNA | Pred unresolved |
| --- | --- | --- | --- | --- |
| rna | 13,989,924 | 13,482,965 (96.38%) | 386,309 (2.76%) | 120,650 (0.86%) |
| gdna | 18,228,677 | 3,331,662 (18.28%) | 14,677,584 (80.52%) | 219,431 (1.20%) |

## Detailed Rigel Pool Calls

| True source | Fragments | Pred mRNA | Pred nRNA | Pred gDNA | Pred unresolved |
| --- | --- | --- | --- | --- | --- |
| rna | 13,989,924 | 12,965,635 (92.68%) | 517,330 (3.70%) | 386,309 (2.76%) | 120,650 (0.86%) |
| gdna | 18,228,677 | 1,385,117 (7.60%) | 1,946,545 (10.68%) | 14,677,584 (80.52%) | 219,431 (1.20%) |

## Largest RNA/gDNA Pool Errors

| Flowcell | Truth | Pred | Detail | ZC | ZS | Count |
| --- | --- | --- | --- | --- | --- | --- |
| H7MFFDSXY | gdna | rna | nrna | ambig_same_strand | unspliced | 1,299,008 |
| H7MFFDSXY | gdna | rna | mrna | ambig_same_strand | unspliced | 1,107,859 |
| H7MFFDSXY | gdna | rna | nrna | ambig_opp_strand | unspliced | 582,324 |
| C6EL5ANXX | rna | gdna | gdna | ambig_same_strand | unspliced | 314,117 |
| H7MFFDSXY | gdna | rna | mrna | ambig_opp_strand | unspliced | 262,978 |
| H7MFFDSXY | gdna | unresolved | unresolved | missing | missing | 147,858 |
| C6EL5ANXX | rna | unresolved | unresolved | missing | missing | 98,133 |
| H7MFFDSXY | gdna | unresolved | unresolved | . | unspliced | 67,776 |
| C6EL5ANXX | rna | gdna | gdna | ambig_opp_strand | unspliced | 48,139 |
| H7MFFDSXY | gdna | rna | nrna | unambig | unspliced | 35,519 |
| H7MFFDSXY | gdna | rna | nrna | multimapper | unspliced | 20,898 |
| C6EL5ANXX | rna | gdna | gdna | unambig | unspliced | 11,292 |
| C6EL5ANXX | rna | unresolved | unresolved | . | spliced_annot | 10,890 |
| C6EL5ANXX | rna | unresolved | unresolved | . | unspliced | 10,387 |
| C6EL5ANXX | rna | gdna | gdna | . | unknown | 9,518 |
| H7MFFDSXY | gdna | rna | mrna | multimapper | unspliced | 7,705 |
| H7MFFDSXY | gdna | rna | nrna | ambig_same_strand | spliced_unannot | 3,773 |
| H7MFFDSXY | gdna | unresolved | unresolved | . | spliced_unannot | 3,704 |
| C6EL5ANXX | rna | gdna | gdna | multimapper | unspliced | 3,243 |
| H7MFFDSXY | gdna | rna | nrna | multimapper | spliced_unannot | 2,694 |

## Artifacts

- Flowcell summary: `results/vcap_rna20m_gdna20m_confusion_2026-05-15/flowcell_summary.tsv`
- RNA/gDNA confusion counts: `results/vcap_rna20m_gdna20m_confusion_2026-05-15/confusion_counts.tsv`
- RNA/gDNA confusion matrix: `results/vcap_rna20m_gdna20m_confusion_2026-05-15/confusion_matrix.tsv`
- Detailed predicted pools: `results/vcap_rna20m_gdna20m_confusion_2026-05-15/detailed_predicted_pool.tsv`
- ZF by flowcell: `results/vcap_rna20m_gdna20m_confusion_2026-05-15/zf_by_flowcell.tsv`
- Error ZC/ZS breakdown: `results/vcap_rna20m_gdna20m_confusion_2026-05-15/error_zc_zs_breakdown.tsv`
- Summary JSON: `results/vcap_rna20m_gdna20m_confusion_2026-05-15/summary.json`
