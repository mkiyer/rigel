# VCaP RNA/gDNA Real Mix Pool Confusion - 2026-05-19

BAM: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/v4_3_with_mm/annotated.bam`

Counting rule: one primary read1 record per fragment; secondary, supplementary, and read2 records are skipped. Truth source is derived from the query-name flowcell ID.

Source mapping:

- `C6EL5ANXX` -> `rna`
- `H7MFFDSXY` -> `gdna`

Fragments counted: 32,218,601

## Headline

Rigel reports 45.29% gDNA and 53.41% RNA across all counted fragments, with 1.31% unresolved. Among assigned RNA/gDNA fragments only, the reported composition is 45.89% gDNA and 54.11% RNA.

Using the flowcell labels as truth, RNA-source fragments are called RNA at 96.37% and gDNA at 2.74%. gDNA-source fragments are called gDNA at 77.94% and RNA at 20.44%.

## Flowcell Summary

| Flowcell | Truth | Fragments | Mean read len | Pred RNA | Pred gDNA | Pred unresolved |
| --- | --- | --- | --- | --- | --- | --- |
| C6EL5ANXX | rna | 13,989,924 | 124.6 | 96.37% | 2.74% | 0.89% |
| H7MFFDSXY | gdna | 18,228,677 | 149.4 | 20.44% | 77.94% | 1.63% |

## RNA versus gDNA Confusion Matrix

| True source | Fragments | Pred RNA | Pred gDNA | Pred unresolved |
| --- | --- | --- | --- | --- |
| rna | 13,989,924 | 13,481,768 (96.37%) | 383,809 (2.74%) | 124,347 (0.89%) |
| gdna | 18,228,677 | 3,725,095 (20.44%) | 14,207,120 (77.94%) | 296,462 (1.63%) |

## Detailed Rigel Pool Calls

| True source | Fragments | Pred mRNA | Pred nRNA | Pred gDNA | Pred unresolved |
| --- | --- | --- | --- | --- | --- |
| rna | 13,989,924 | 12,921,620 (92.36%) | 560,148 (4.00%) | 383,809 (2.74%) | 124,347 (0.89%) |
| gdna | 18,228,677 | 1,346,710 (7.39%) | 2,378,385 (13.05%) | 14,207,120 (77.94%) | 296,462 (1.63%) |

## Largest RNA/gDNA Pool Errors

| Flowcell | Truth | Pred | Detail | ZC | ZS | Count |
| --- | --- | --- | --- | --- | --- | --- |
| H7MFFDSXY | gdna | rna | nrna | ambig_same_strand | unspliced | 1,661,933 |
| H7MFFDSXY | gdna | rna | mrna | ambig_same_strand | unspliced | 1,048,449 |
| H7MFFDSXY | gdna | rna | nrna | ambig_opp_strand | unspliced | 598,628 |
| C6EL5ANXX | rna | gdna | gdna | ambig_same_strand | unspliced | 311,221 |
| H7MFFDSXY | gdna | rna | mrna | ambig_opp_strand | unspliced | 260,025 |
| H7MFFDSXY | gdna | unresolved | unresolved | missing | missing | 147,858 |
| H7MFFDSXY | gdna | unresolved | unresolved | . | unspliced | 145,027 |
| C6EL5ANXX | rna | unresolved | unresolved | missing | missing | 98,133 |
| C6EL5ANXX | rna | gdna | gdna | ambig_opp_strand | unspliced | 47,211 |
| H7MFFDSXY | gdna | rna | nrna | ambig_same_strand | spliced_implicit | 41,217 |
| H7MFFDSXY | gdna | rna | nrna | multimapper | unspliced | 25,463 |
| H7MFFDSXY | gdna | rna | nrna | unambig | unspliced | 25,266 |
| C6EL5ANXX | rna | gdna | gdna | unambig | unspliced | 13,330 |
| C6EL5ANXX | rna | unresolved | unresolved | . | spliced_annot | 11,966 |
| C6EL5ANXX | rna | unresolved | unresolved | . | unspliced | 11,703 |
| H7MFFDSXY | gdna | rna | nrna | multimapper | spliced_annot | 10,525 |
| C6EL5ANXX | rna | gdna | gdna | . | unknown | 9,491 |
| H7MFFDSXY | gdna | rna | mrna | ambig_same_strand | spliced_annot | 9,159 |
| H7MFFDSXY | gdna | rna | mrna | multimapper | unspliced | 8,364 |
| H7MFFDSXY | gdna | rna | mrna | ambig_same_strand | spliced_implicit | 7,302 |

## Artifacts

- Flowcell summary: `results/vcap_rna20m_gdna20m_v4_3_confusion_2026-05-19/flowcell_summary.tsv`
- RNA/gDNA confusion counts: `results/vcap_rna20m_gdna20m_v4_3_confusion_2026-05-19/confusion_counts.tsv`
- RNA/gDNA confusion matrix: `results/vcap_rna20m_gdna20m_v4_3_confusion_2026-05-19/confusion_matrix.tsv`
- Detailed predicted pools: `results/vcap_rna20m_gdna20m_v4_3_confusion_2026-05-19/detailed_predicted_pool.tsv`
- ZF by flowcell: `results/vcap_rna20m_gdna20m_v4_3_confusion_2026-05-19/zf_by_flowcell.tsv`
- Error ZC/ZS breakdown: `results/vcap_rna20m_gdna20m_v4_3_confusion_2026-05-19/error_zc_zs_breakdown.tsv`
- Summary JSON: `results/vcap_rna20m_gdna20m_v4_3_confusion_2026-05-19/summary.json`
